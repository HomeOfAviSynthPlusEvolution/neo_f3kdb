// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <int kDitherAlgo>
int postprocess_scalar_pixel(
    const neo_f3kdb::core::KernelConfig& config,
    int pixel,
    span2d::ReadOnlyRestrictRow<std::int16_t> grain_row,
    int row,
    int col
) {
    pixel += static_cast<int>(grain_row[static_cast<std::size_t>(col)]);

    if constexpr (kDitherAlgo == DA_HIGH_NO_DITHERING) {
        return std::clamp(pixel, config.pixel_min, config.pixel_max) >> (16 - config.output_depth);
    } else if constexpr (kDitherAlgo == DA_HIGH_ORDERED_DITHERING) {
        pixel += neo_f3kdb::core::pixel_proc_detail::ordered::THRESHOLD_MAP[row & 15][col & 15] >> (config.output_depth - 8);
        return std::clamp(pixel, config.pixel_min, config.pixel_max) >> (16 - config.output_depth);
    } else {
        static_assert(kDitherAlgo == DA_16BIT_INTERLEAVED);
        return std::clamp(pixel, config.pixel_min, config.pixel_max);
    }
}

inline int postprocess_floyd_steinberg_pixel(
    const neo_f3kdb::core::KernelConfig& config,
    neo_f3kdb::core::dither::FloydSteinbergDither& fs_dither,
    int pixel,
    span2d::ReadOnlyRestrictRow<std::int16_t> grain_row,
    int col
) {
    pixel += static_cast<int>(grain_row[static_cast<std::size_t>(col)]);
    pixel = fs_dither.dither(pixel);
    return neo_f3kdb::core::pixel_proc_common::downsample(
        pixel,
        config.pixel_min,
        config.pixel_max,
        config.output_depth
    );
}

template <int kDitherAlgo, class D32, class V32>
V32 apply_dither_and_grain(
    const neo_f3kdb::core::KernelConfig& config,
    D32 d32,
    V32 pixel,
    span2d::ReadOnlyRestrictRow<std::int16_t> grain_row,
    int row,
    int col,
    std::size_t lanes
) {
    const hn::Rebind<std::int16_t, D32> d16;
    pixel = hn::Add(pixel, hn::PromoteTo(d32, hn::LoadU(d16, grain_row.data() + col)));

    if constexpr (kDitherAlgo == DA_HIGH_ORDERED_DITHERING) {
        alignas(64) std::int32_t dither[hn::MaxLanes(d32)];
        for (std::size_t lane = 0; lane < lanes; ++lane) {
            dither[lane] = neo_f3kdb::core::pixel_proc_detail::ordered::THRESHOLD_MAP[row & 15][(col + static_cast<int>(lane)) & 15] >> (config.output_depth - 8);
        }
        pixel = hn::Add(pixel, hn::LoadU(d32, dither));
    }

    pixel = hn::Min(hn::Max(pixel, hn::Set(d32, config.pixel_min)), hn::Set(d32, config.pixel_max));
    if constexpr (kDitherAlgo != DA_16BIT_INTERLEAVED) {
        pixel = hn::ShiftRightSame(pixel, 16 - config.output_depth);
    }
    return pixel;
}

template <class DstType>
HWY_INLINE void process_stripe_fs_wavefront_hwy(
    const neo_f3kdb::core::KernelConfig& config,
    const std::uint16_t* const* stripe_src_rows,
    DstType* const* stripe_dst_rows,
    int width,
    int height,
    int y0,
    int stripe_height,
    std::vector<std::uint16_t>& cur_bridge,
    std::vector<std::uint16_t>& next_bridge
) {
    const hn::ScalableTag<std::uint16_t> d;
    const std::size_t N = hn::Lanes(d);

    const int shift = 16 - config.output_depth;
    const std::uint16_t err_mask = static_cast<std::uint16_t>((1 << shift) - 1);

    const auto v_mask = hn::Set(d, err_mask);
    const auto v_min = hn::Set(d, static_cast<std::uint16_t>(config.pixel_min));
    const auto v_max = hn::Set(d, static_cast<std::uint16_t>(config.pixel_max));

    const auto v_w_left = hn::Set(d, 7);
    const auto v_w_tr   = hn::Set(d, 3);
    const auto v_w_top  = hn::Set(d, 5);
    const auto v_w_tl   = hn::Set(d, 1);

    alignas(64) std::uint16_t src_arr[hn::MaxLanes(d)];
    alignas(64) std::uint16_t dst_arr[hn::MaxLanes(d)];
    alignas(64) std::uint16_t err_arr[hn::MaxLanes(d)];

    auto V_prev1 = hn::Zero(d);
    auto V_prev2 = hn::Zero(d);
    auto V_prev3 = hn::Zero(d);

    const int max_t = width + 2 * static_cast<int>(N) - 2;

    if (stripe_height == static_cast<int>(N)) {
        const int t_steady_start = 2 * (static_cast<int>(N) - 1);
        const int t_steady_end = width;

        // Phase 1: Ramp-up (t < t_steady_start)
        for (int t = 0; t < t_steady_start && t < max_t; ++t) {
            auto v_left = hn::ShiftRight<4>(hn::Mul(V_prev1, v_w_left));
            auto v_tr   = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev1), v_w_tr));
            auto v_top  = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev2), v_w_top));
            auto v_tl   = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev3), v_w_tl));
            auto v_err  = hn::Add(hn::Add(v_left, v_tr), hn::Add(v_top, v_tl));

            if (t >= 0 && t < width && y0 > 0) {
                v_err = hn::Add(v_err, hn::InsertLane(hn::Zero(d), 0, cur_bridge[static_cast<std::size_t>(t)]));
            }

            for (std::size_t k = 0; k < N; ++k) {
                int xk = t - 2 * static_cast<int>(k);
                src_arr[k] = (xk >= 0 && xk < width) ? stripe_src_rows[k][xk] : 0;
            }

            auto v_src = hn::LoadU(d, src_arr);
            auto v_pix = hn::SaturatedAdd(v_src, v_err);
            auto v_new_err = hn::And(v_pix, v_mask);
            auto v_clamped = hn::Min(hn::Max(v_pix, v_min), v_max);
            auto v_out = hn::ShiftRightSame(v_clamped, shift);

            hn::StoreU(v_out, d, dst_arr);
            hn::StoreU(v_new_err, d, err_arr);

            for (std::size_t k = 0; k < N; ++k) {
                int xk = t - 2 * static_cast<int>(k);
                if (xk >= 0 && xk < width) {
                    stripe_dst_rows[k][xk] = static_cast<DstType>(dst_arr[k]);
                } else {
                    err_arr[k] = 0;
                }
            }

            const int last_k = static_cast<int>(N) - 1;
            const int x_last = t - 2 * last_k;
            if (x_last >= 0 && x_last < width && y0 + static_cast<int>(N) < height) {
                const std::uint16_t e = err_arr[last_k];
                if (x_last - 1 >= 0)    next_bridge[static_cast<std::size_t>(x_last - 1)] += (e * 3) >> 4;
                if (x_last < width)     next_bridge[static_cast<std::size_t>(x_last)]     += (e * 5) >> 4;
                if (x_last + 1 < width) next_bridge[static_cast<std::size_t>(x_last + 1)] += (e * 1) >> 4;
            }

            V_prev3 = V_prev2;
            V_prev2 = V_prev1;
            V_prev1 = hn::LoadU(d, err_arr);
        }

        // Phase 2: Steady State (all N lanes active, 100% in-register, pointer streaming)
        const std::uint16_t* ptr_src[hn::MaxLanes(d)];
        DstType* ptr_dst[hn::MaxLanes(d)];
        for (std::size_t k = 0; k < N; ++k) {
            ptr_src[k] = stripe_src_rows[k] + (t_steady_start - 2 * static_cast<int>(k));
            ptr_dst[k] = stripe_dst_rows[k] + (t_steady_start - 2 * static_cast<int>(k));
        }

        for (int t = t_steady_start; t < t_steady_end; ++t) {
            auto v_left = hn::ShiftRight<4>(hn::Mul(V_prev1, v_w_left));
            auto v_tr   = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev1), v_w_tr));
            auto v_top  = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev2), v_w_top));
            auto v_tl   = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev3), v_w_tl));
            auto v_err  = hn::Add(hn::Add(v_left, v_tr), hn::Add(v_top, v_tl));

            if (y0 > 0) {
                v_err = hn::Add(v_err, hn::InsertLane(hn::Zero(d), 0, cur_bridge[static_cast<std::size_t>(t)]));
            }

            for (std::size_t k = 0; k < N; ++k) {
                src_arr[k] = *ptr_src[k]++;
            }

            auto v_src = hn::LoadU(d, src_arr);
            auto v_pix = hn::SaturatedAdd(v_src, v_err);
            auto v_new_err = hn::And(v_pix, v_mask);
            auto v_clamped = hn::Min(hn::Max(v_pix, v_min), v_max);
            auto v_out = hn::ShiftRightSame(v_clamped, shift);

            hn::StoreU(v_out, d, dst_arr);

            for (std::size_t k = 0; k < N; ++k) {
                *ptr_dst[k]++ = static_cast<DstType>(dst_arr[k]);
            }

            if (y0 + static_cast<int>(N) < height) {
                const int last_k = static_cast<int>(N) - 1;
                const int x_last = t - 2 * last_k;
                const std::uint16_t e = hn::ExtractLane(v_new_err, last_k);
                if (x_last - 1 >= 0)    next_bridge[static_cast<std::size_t>(x_last - 1)] += (e * 3) >> 4;
                if (x_last < width)     next_bridge[static_cast<std::size_t>(x_last)]     += (e * 5) >> 4;
                if (x_last + 1 < width) next_bridge[static_cast<std::size_t>(x_last + 1)] += (e * 1) >> 4;
            }

            V_prev3 = V_prev2;
            V_prev2 = V_prev1;
            V_prev1 = v_new_err; // 100% in-register feedback!
        }

        // Phase 3: Ramp-down (t >= t_steady_end)
        for (int t = t_steady_end; t < max_t; ++t) {
            auto v_left = hn::ShiftRight<4>(hn::Mul(V_prev1, v_w_left));
            auto v_tr   = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev1), v_w_tr));
            auto v_top  = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev2), v_w_top));
            auto v_tl   = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev3), v_w_tl));
            auto v_err  = hn::Add(hn::Add(v_left, v_tr), hn::Add(v_top, v_tl));

            for (std::size_t k = 0; k < N; ++k) {
                int xk = t - 2 * static_cast<int>(k);
                src_arr[k] = (xk >= 0 && xk < width) ? stripe_src_rows[k][xk] : 0;
            }

            auto v_src = hn::LoadU(d, src_arr);
            auto v_pix = hn::SaturatedAdd(v_src, v_err);
            auto v_new_err = hn::And(v_pix, v_mask);
            auto v_clamped = hn::Min(hn::Max(v_pix, v_min), v_max);
            auto v_out = hn::ShiftRightSame(v_clamped, shift);

            hn::StoreU(v_out, d, dst_arr);
            hn::StoreU(v_new_err, d, err_arr);

            for (std::size_t k = 0; k < N; ++k) {
                int xk = t - 2 * static_cast<int>(k);
                if (xk >= 0 && xk < width) {
                    stripe_dst_rows[k][xk] = static_cast<DstType>(dst_arr[k]);
                } else {
                    err_arr[k] = 0;
                }
            }

            const int last_k = static_cast<int>(N) - 1;
            const int x_last = t - 2 * last_k;
            if (x_last >= 0 && x_last < width && y0 + static_cast<int>(N) < height) {
                const std::uint16_t e = err_arr[last_k];
                if (x_last - 1 >= 0)    next_bridge[static_cast<std::size_t>(x_last - 1)] += (e * 3) >> 4;
                if (x_last < width)     next_bridge[static_cast<std::size_t>(x_last)]     += (e * 5) >> 4;
                if (x_last + 1 < width) next_bridge[static_cast<std::size_t>(x_last + 1)] += (e * 1) >> 4;
            }

            V_prev3 = V_prev2;
            V_prev2 = V_prev1;
            V_prev1 = hn::LoadU(d, err_arr);
        }
    } else {
        // Tail stripe fallback (when stripe_height < N)
        for (int t = 0; t < max_t; ++t) {
            auto v_left = hn::ShiftRight<4>(hn::Mul(V_prev1, v_w_left));
            auto v_tr   = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev1), v_w_tr));
            auto v_top  = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev2), v_w_top));
            auto v_tl   = hn::ShiftRight<4>(hn::Mul(hn::Slide1Up(d, V_prev3), v_w_tl));
            auto v_err  = hn::Add(hn::Add(v_left, v_tr), hn::Add(v_top, v_tl));

            if (t >= 0 && t < width && y0 > 0) {
                v_err = hn::Add(v_err, hn::InsertLane(hn::Zero(d), 0, cur_bridge[static_cast<std::size_t>(t)]));
            }

            for (int k = 0; k < stripe_height; ++k) {
                int xk = t - 2 * k;
                src_arr[k] = (xk >= 0 && xk < width) ? stripe_src_rows[k][xk] : 0;
            }
            for (int k = stripe_height; k < static_cast<int>(N); ++k) {
                src_arr[k] = 0;
            }

            auto v_src = hn::LoadU(d, src_arr);
            auto v_pix = hn::SaturatedAdd(v_src, v_err);
            auto v_new_err = hn::And(v_pix, v_mask);
            auto v_clamped = hn::Min(hn::Max(v_pix, v_min), v_max);
            auto v_out = hn::ShiftRightSame(v_clamped, shift);

            hn::StoreU(v_out, d, dst_arr);
            hn::StoreU(v_new_err, d, err_arr);

            for (int k = 0; k < stripe_height; ++k) {
                int xk = t - 2 * k;
                if (xk >= 0 && xk < width) {
                    stripe_dst_rows[k][xk] = static_cast<DstType>(dst_arr[k]);
                } else {
                    err_arr[k] = 0;
                }
            }
            for (int k = stripe_height; k < static_cast<int>(N); ++k) {
                err_arr[k] = 0;
            }

            const int last_k = stripe_height - 1;
            const int x_last = t - 2 * last_k;
            if (x_last >= 0 && x_last < width && y0 + stripe_height < height) {
                const std::uint16_t e = err_arr[last_k];
                if (x_last - 1 >= 0)    next_bridge[static_cast<std::size_t>(x_last - 1)] += (e * 3) >> 4;
                if (x_last < width)     next_bridge[static_cast<std::size_t>(x_last)]     += (e * 5) >> 4;
                if (x_last + 1 < width) next_bridge[static_cast<std::size_t>(x_last + 1)] += (e * 1) >> 4;
            }

            V_prev3 = V_prev2;
            V_prev2 = V_prev1;
            V_prev1 = hn::LoadU(d, err_arr);
        }
    }
}
