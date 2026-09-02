// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <int kSampleMode, class PixelIn>
HWY_NOINLINE neo_f3kdb::core::PlaneOffsetCache* build_offset_cache(
    const process_plane_params& params,
    span2d::ReadOnlyRestrictPlane<PixelIn> src_plane
) {
    static_assert(kSampleMode >= 1 && kSampleMode <= 6);

    const int width = params.plane_width();
    const int height = params.plane_height();
    auto* cache = new neo_f3kdb::core::PlaneOffsetCache();
    cache->pitch = static_cast<int>(src_plane.stride_bytes());
    cache->width = width;
    cache->height = height;
    cache->sample_mode = kSampleMode;
    cache->input_depth = params.config.input_depth;

    const std::size_t total_pixels = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    cache->off1.resize(total_pixels);
    cache->off2.resize(total_pixels);
    cache->off3.resize(total_pixels);
    cache->off4.resize(total_pixels);

    const int max_y = height - 1;
    const int max_x = width - 1;
    const int width_subsamp = params.config.width_subsampling;
    const int height_subsamp = params.config.height_subsampling;
    const int byte_step = static_cast<int>(sizeof(PixelIn));

    for (int row = 0; row < height; ++row) {
        const auto* dinfo = params.dither_info_plane().row_ptr(row);
        const auto* row_base = reinterpret_cast<const unsigned char*>(src_plane.row_ptr(row));
        for (int col = 0; col < width; ++col) {
            const std::size_t idx = static_cast<std::size_t>(row) * width + col;
            const auto* current = row_base + col * byte_step;
            const auto& info = dinfo[col];
            const int ref_y1 = info.ref1 >> height_subsamp;
            const int ref_y2 = info.ref2 >> height_subsamp;
            const int ref_x1 = info.ref1 >> width_subsamp;
            const int ref_x2 = info.ref2 >> width_subsamp;

            int y1 = row, x1 = col, y2 = row, x2 = col, y3 = row, x3 = col, y4 = row, x4 = col;
            if constexpr (kSampleMode == 1) {
                y1 = std::clamp(row + ref_y1, 0, max_y); x1 = col;
                y2 = std::clamp(row - ref_y1, 0, max_y); x2 = col;
            } else if constexpr (kSampleMode == 2) {
                y1 = std::clamp(row + ref_y2, 0, max_y); x1 = std::clamp(col + ref_x1, 0, max_x);
                y2 = std::clamp(row - ref_y1, 0, max_y); x2 = std::clamp(col + ref_x2, 0, max_x);
                y3 = std::clamp(row - ref_y2, 0, max_y); x3 = std::clamp(col - ref_x1, 0, max_x);
                y4 = std::clamp(row + ref_y1, 0, max_y); x4 = std::clamp(col - ref_x2, 0, max_x);
            } else if constexpr (kSampleMode == 3) {
                y1 = row; x1 = std::clamp(col + ref_x1, 0, max_x);
                y2 = row; x2 = std::clamp(col - ref_x1, 0, max_x);
            } else if constexpr (kSampleMode == 4 || kSampleMode == 5 || kSampleMode == 6) {
                y1 = std::clamp(row + ref_y1, 0, max_y); x1 = col;
                y2 = std::clamp(row - ref_y1, 0, max_y); x2 = col;
                y3 = row; x3 = std::clamp(col + ref_x1, 0, max_x);
                y4 = row; x4 = std::clamp(col - ref_x1, 0, max_x);
            }

            const auto* p1 = reinterpret_cast<const unsigned char*>(src_plane.row_ptr(y1)) + x1 * byte_step;
            const auto* p2 = reinterpret_cast<const unsigned char*>(src_plane.row_ptr(y2)) + x2 * byte_step;
            const auto* p3 = reinterpret_cast<const unsigned char*>(src_plane.row_ptr(y3)) + x3 * byte_step;
            const auto* p4 = reinterpret_cast<const unsigned char*>(src_plane.row_ptr(y4)) + x4 * byte_step;

            cache->off1[idx] = static_cast<std::int32_t>(p1 - current);
            cache->off2[idx] = static_cast<std::int32_t>(p2 - current);
            cache->off3[idx] = static_cast<std::int32_t>(p3 - current);
            cache->off4[idx] = static_cast<std::int32_t>(p4 - current);
        }
    }

    return cache;
}

template <int kSampleMode, class PixelIn>
HWY_NOINLINE const neo_f3kdb::core::PlaneOffsetCache* prepare_offset_cache(
    const process_plane_params& params,
    span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
    process_plane_context* context
) {
    static_assert(kSampleMode >= 1 && kSampleMode <= 6);

    if (!context) {
        return nullptr;
    }

    const int pitch = static_cast<int>(src_plane.stride_bytes());
    const int width = params.plane_width();
    const int height = params.plane_height();
    const int depth = params.config.input_depth;

    // Fast path: lock-free check with acquire semantics
    void* existing_ptr = context->data.load(std::memory_order_acquire);
    if (existing_ptr) {
        const auto* cache = static_cast<const neo_f3kdb::core::PlaneOffsetCache*>(existing_ptr);
        if (cache && cache->matches(pitch, width, height, kSampleMode, depth)) {
            return cache;
        }
    }

    // Slow path: lock mutex and perform double-checked initialization
    std::lock_guard<std::mutex> lock(context->mutex);
    existing_ptr = context->data.load(std::memory_order_relaxed);
    if (existing_ptr) {
        const auto* cache = static_cast<const neo_f3kdb::core::PlaneOffsetCache*>(existing_ptr);
        if (cache && cache->matches(pitch, width, height, kSampleMode, depth)) {
            return cache;
        }
        if (context->destroy) {
            context->destroy(existing_ptr);
        }
    }

    auto* cache = build_offset_cache<kSampleMode>(params, src_plane);
    context->destroy = neo_f3kdb::core::destroy_offset_cache;
    context->data.store(cache, std::memory_order_release);
    return cache;
}

template <class DTag, class PixelIn, class VTag>
HWY_INLINE void gather_cached_reference_vectors(
    DTag dtag,
    const PixelIn* current,
    const std::int32_t* off1,
    const std::int32_t* off2,
    const std::int32_t* off3,
    const std::int32_t* off4,
    int upshift,
    std::size_t lanes,
    VTag& r1_v,
    VTag& r2_v,
    VTag& r3_v,
    VTag& r4_v,
    bool is_chroma = false
) {
    using T = hn::TFromD<DTag>;
    alignas(64) T ref1[hn::MaxLanes(dtag)];
    alignas(64) T ref2[hn::MaxLanes(dtag)];
    alignas(64) T ref3[hn::MaxLanes(dtag)];
    alignas(64) T ref4[hn::MaxLanes(dtag)];

    const auto* current_bytes = reinterpret_cast<const unsigned char*>(current);
    for (std::size_t lane = 0; lane < lanes; ++lane) {
        const auto* lane_bytes = current_bytes + lane * sizeof(PixelIn);
        const PixelIn p1 = *reinterpret_cast<const PixelIn*>(lane_bytes + off1[lane]);
        const PixelIn p2 = *reinterpret_cast<const PixelIn*>(lane_bytes + off2[lane]);
        const PixelIn p3 = *reinterpret_cast<const PixelIn*>(lane_bytes + off3[lane]);
        const PixelIn p4 = *reinterpret_cast<const PixelIn*>(lane_bytes + off4[lane]);
        ref1[lane] = static_cast<T>(neo_f3kdb::core::sample_modes::upsample(p1, upshift == 0 ? 32 : (16 - upshift), is_chroma));
        ref2[lane] = static_cast<T>(neo_f3kdb::core::sample_modes::upsample(p2, upshift == 0 ? 32 : (16 - upshift), is_chroma));
        ref3[lane] = static_cast<T>(neo_f3kdb::core::sample_modes::upsample(p3, upshift == 0 ? 32 : (16 - upshift), is_chroma));
        ref4[lane] = static_cast<T>(neo_f3kdb::core::sample_modes::upsample(p4, upshift == 0 ? 32 : (16 - upshift), is_chroma));
    }

    r1_v = hn::LoadU(dtag, ref1);
    r2_v = hn::LoadU(dtag, ref2);
    r3_v = hn::LoadU(dtag, ref3);
    r4_v = hn::LoadU(dtag, ref4);
}

template <int kSampleMode, class DTag, class PixelIn, class VTag>
HWY_INLINE void gather_uncached_reference_vectors(
    DTag dtag,
    const process_plane_params& params,
    span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
    int row,
    int col,
    std::size_t lanes,
    VTag& r1_v,
    VTag& r2_v,
    VTag& r3_v,
    VTag& r4_v
) {
    using T = hn::TFromD<DTag>;
    static_assert(kSampleMode >= 1 && kSampleMode <= 7);

    alignas(64) T ref1[hn::MaxLanes(dtag)];
    alignas(64) T ref2[hn::MaxLanes(dtag)];
    alignas(64) T ref3[hn::MaxLanes(dtag)];
    alignas(64) T ref4[hn::MaxLanes(dtag)];

    const auto* dinfo = params.dither_info_plane().row_ptr(row) + col;
    const int width_subsamp = params.config.width_subsampling;
    const int height_subsamp = params.config.height_subsampling;
    const int max_y = params.plane_height() - 1;
    const int max_x = params.plane_width() - 1;

    if constexpr (kSampleMode == 2) {
        for (std::size_t lane = 0; lane < lanes; ++lane) {
            const int pixel_col = col + static_cast<int>(lane);
            const auto& info = dinfo[lane];
            const int ref_y1 = info.ref1 >> height_subsamp;
            const int ref_y2 = info.ref2 >> height_subsamp;
            const int ref_x1 = info.ref1 >> width_subsamp;
            const int ref_x2 = info.ref2 >> width_subsamp;

            const bool is_chroma = params.config.plane_index > 0;
            const int y1 = std::clamp(row + ref_y2, 0, max_y);
            const int x1 = std::clamp(pixel_col + ref_x1, 0, max_x);
            ref1[lane] = static_cast<T>(neo_f3kdb::core::sample_modes::upsample(src_plane(y1, x1), params.config.input_depth, is_chroma));

            const int y2 = std::clamp(row - ref_y1, 0, max_y);
            const int x2 = std::clamp(pixel_col + ref_x2, 0, max_x);
            ref2[lane] = static_cast<T>(neo_f3kdb::core::sample_modes::upsample(src_plane(y2, x2), params.config.input_depth, is_chroma));

            const int y3 = std::clamp(row - ref_y2, 0, max_y);
            const int x3 = std::clamp(pixel_col - ref_x1, 0, max_x);
            ref3[lane] = static_cast<T>(neo_f3kdb::core::sample_modes::upsample(src_plane(y3, x3), params.config.input_depth, is_chroma));

            const int y4 = std::clamp(row + ref_y1, 0, max_y);
            const int x4 = std::clamp(pixel_col - ref_x2, 0, max_x);
            ref4[lane] = static_cast<T>(neo_f3kdb::core::sample_modes::upsample(src_plane(y4, x4), params.config.input_depth, is_chroma));
        }
    } else {
        const neo_f3kdb::core::sample_modes::TypedSampleOps<PixelIn> ops{params, src_plane};
        for (std::size_t lane = 0; lane < lanes; ++lane) {
            const int pixel_col = col + static_cast<int>(lane);
            int r1 = 0, r2 = 0, r3 = 0, r4 = 0;
            neo_f3kdb::core::sample_modes::load_reference_samples_only_ref<kSampleMode>(
                params,
                ops,
                row,
                pixel_col,
                r1, r2, r3, r4
            );
            ref1[lane] = static_cast<T>(r1);
            ref2[lane] = static_cast<T>(r2);
            ref3[lane] = static_cast<T>(r3);
            ref4[lane] = static_cast<T>(r4);
        }
    }

    r1_v = hn::LoadU(dtag, ref1);
    r2_v = hn::LoadU(dtag, ref2);
    r3_v = hn::LoadU(dtag, ref3);
    r4_v = hn::LoadU(dtag, ref4);
}

template <class DF, class PixelIn>
HWY_INLINE hn::Vec<DF> gather_gradient_angle_vector(
    DF df,
    const process_plane_params& params,
    span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
    const int* y_coords,
    const int* x_coords,
    std::size_t lanes,
    float scaled_eps_gx,
    int read_distance = 20
) {
    const int max_y = params.plane_height() - 1;
    const int max_x = params.plane_width() - 1;

    alignas(64) float p00_b[hn::MaxLanes(df)];
    alignas(64) float p10_b[hn::MaxLanes(df)];
    alignas(64) float p20_b[hn::MaxLanes(df)];
    alignas(64) float p01_b[hn::MaxLanes(df)];
    alignas(64) float p21_b[hn::MaxLanes(df)];
    alignas(64) float p02_b[hn::MaxLanes(df)];
    alignas(64) float p12_b[hn::MaxLanes(df)];
    alignas(64) float p22_b[hn::MaxLanes(df)];

    for (std::size_t i = 0; i < lanes; ++i) {
        const int cy = y_coords[i];
        const int cx = x_coords[i];

        const int ym = std::clamp(cy - read_distance, 0, max_y);
        const int y0 = std::clamp(cy, 0, max_y);
        const int yp = std::clamp(cy + read_distance, 0, max_y);

        const int xm = std::clamp(cx - read_distance, 0, max_x);
        const int x0 = std::clamp(cx, 0, max_x);
        const int xp = std::clamp(cx + read_distance, 0, max_x);

        const auto* row_m = src_plane.row_ptr(ym);
        const auto* row_0 = src_plane.row_ptr(y0);
        const auto* row_p = src_plane.row_ptr(yp);

        p00_b[i] = static_cast<float>(neo_f3kdb::core::sample_modes::upsample(row_m[xm], params.config.input_depth));
        p10_b[i] = static_cast<float>(neo_f3kdb::core::sample_modes::upsample(row_m[x0], params.config.input_depth));
        p20_b[i] = static_cast<float>(neo_f3kdb::core::sample_modes::upsample(row_m[xp], params.config.input_depth));

        p01_b[i] = static_cast<float>(neo_f3kdb::core::sample_modes::upsample(row_0[xm], params.config.input_depth));
        p21_b[i] = static_cast<float>(neo_f3kdb::core::sample_modes::upsample(row_0[xp], params.config.input_depth));

        p02_b[i] = static_cast<float>(neo_f3kdb::core::sample_modes::upsample(row_p[xm], params.config.input_depth));
        p12_b[i] = static_cast<float>(neo_f3kdb::core::sample_modes::upsample(row_p[x0], params.config.input_depth));
        p22_b[i] = static_cast<float>(neo_f3kdb::core::sample_modes::upsample(row_p[xp], params.config.input_depth));
    }

    const auto p00 = hn::LoadU(df, p00_b);
    const auto p10 = hn::LoadU(df, p10_b);
    const auto p20 = hn::LoadU(df, p20_b);
    const auto p01 = hn::LoadU(df, p01_b);
    const auto p21 = hn::LoadU(df, p21_b);
    const auto p02 = hn::LoadU(df, p02_b);
    const auto p12 = hn::LoadU(df, p12_b);
    const auto p22 = hn::LoadU(df, p22_b);

    const auto two = hn::Set(df, 2.0f);
    const auto inv_pi = hn::Set(df, 1.0f / 3.14159265358979323846f);
    const auto half = hn::Set(df, 0.5f);

    const auto gx = hn::Sub(hn::Add(hn::Add(p20, hn::Mul(two, p21)), p22), hn::Add(hn::Add(p00, hn::Mul(two, p01)), p02));
    const auto gy = hn::Sub(hn::Add(hn::Add(p00, hn::Mul(two, p10)), p20), hn::Add(hn::Add(p02, hn::Mul(two, p12)), p22));

    const auto is_zero = hn::Lt(hn::Abs(gx), hn::Set(df, scaled_eps_gx));
    const auto safe_gx = hn::IfThenElse(is_zero, hn::Set(df, 1.0f), gx);
    const auto atan_val = hn::Atan(df, hn::Div(gy, safe_gx));
    const auto angle = hn::Add(hn::Mul(atan_val, inv_pi), half);

    return hn::IfThenElse(is_zero, hn::Set(df, 1.0f), angle);
}
