// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <int kSampleMode, class PixelIn>
HWY_NOINLINE neo_f3kdb::core::PlaneOffsetCache* build_offset_cache(
    const process_plane_params& params,
    span2d::ReadOnlyRestrictPlane<PixelIn> src_plane
) {
    static_assert(kSampleMode >= 1 && kSampleMode <= 5);

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
            } else if constexpr (kSampleMode == 4 || kSampleMode == 5) {
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
    static_assert(kSampleMode >= 1 && kSampleMode <= 5);

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
    VTag& r4_v
) {
    using T = hn::TFromD<DTag>;
    alignas(64) T ref1[hn::MaxLanes(dtag)];
    alignas(64) T ref2[hn::MaxLanes(dtag)];
    alignas(64) T ref3[hn::MaxLanes(dtag)];
    alignas(64) T ref4[hn::MaxLanes(dtag)];

    const auto* current_bytes = reinterpret_cast<const unsigned char*>(current);
    for (std::size_t lane = 0; lane < lanes; ++lane) {
        const auto* lane_bytes = current_bytes + lane * sizeof(PixelIn);
        ref1[lane] = static_cast<T>(static_cast<int>(*reinterpret_cast<const PixelIn*>(lane_bytes + off1[lane])) << upshift);
        ref2[lane] = static_cast<T>(static_cast<int>(*reinterpret_cast<const PixelIn*>(lane_bytes + off2[lane])) << upshift);
        ref3[lane] = static_cast<T>(static_cast<int>(*reinterpret_cast<const PixelIn*>(lane_bytes + off3[lane])) << upshift);
        ref4[lane] = static_cast<T>(static_cast<int>(*reinterpret_cast<const PixelIn*>(lane_bytes + off4[lane])) << upshift);
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
    static_assert(kSampleMode >= 1 && kSampleMode <= 5);

    const int upshift = 16 - params.config.input_depth;
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

            const int y1 = std::clamp(row + ref_y2, 0, max_y);
            const int x1 = std::clamp(pixel_col + ref_x1, 0, max_x);
            ref1[lane] = static_cast<T>(static_cast<int>(src_plane(y1, x1)) << upshift);

            const int y2 = std::clamp(row - ref_y1, 0, max_y);
            const int x2 = std::clamp(pixel_col + ref_x2, 0, max_x);
            ref2[lane] = static_cast<T>(static_cast<int>(src_plane(y2, x2)) << upshift);

            const int y3 = std::clamp(row - ref_y2, 0, max_y);
            const int x3 = std::clamp(pixel_col - ref_x1, 0, max_x);
            ref3[lane] = static_cast<T>(static_cast<int>(src_plane(y3, x3)) << upshift);

            const int y4 = std::clamp(row + ref_y1, 0, max_y);
            const int x4 = std::clamp(pixel_col - ref_x2, 0, max_x);
            ref4[lane] = static_cast<T>(static_cast<int>(src_plane(y4, x4)) << upshift);
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
