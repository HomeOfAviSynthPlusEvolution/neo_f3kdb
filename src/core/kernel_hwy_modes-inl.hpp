// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <int kSampleMode, bool kBlurFirst, class D16, class V16>
V16 process_reference_samples_u16(
    D16 d16,
    V16 src,
    V16 r1,
    V16 r2,
    V16 r3,
    V16 r4,
    std::uint16_t threshold
) {
    static_assert(kSampleMode >= 1 && kSampleMode <= 4);

    const auto threshold_v = hn::Set(d16, threshold);

    if constexpr (kSampleMode == 1 || kSampleMode == 3) {
        const auto avg = hn::AverageRound(r1, r2);
        const auto use_src = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, r1), threshold_v), hn::Ge(hn::AbsDiff(src, r2), threshold_v));
        return hn::IfThenElse(use_src, src, avg);
    } else if constexpr (kSampleMode == 2) {
        auto avg1 = hn::AverageRound(r1, r2);
        const auto avg2 = hn::AverageRound(r3, r4);
        avg1 = hn::SaturatedSub(avg1, hn::Set(d16, static_cast<std::uint16_t>(1)));
        const auto avg = hn::AverageRound(avg1, avg2);
        const auto use_src = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg, src), threshold_v)
            : hn::Or(
                hn::Or(hn::Ge(hn::AbsDiff(r1, src), threshold_v), hn::Ge(hn::AbsDiff(r2, src), threshold_v)),
                hn::Or(hn::Ge(hn::AbsDiff(r3, src), threshold_v), hn::Ge(hn::AbsDiff(r4, src), threshold_v))
            );
        return hn::IfThenElse(use_src, src, avg);
    } else {
        static_assert(kSampleMode == 4);
        const auto avg_v = hn::AverageRound(r1, r2);
        const auto avg_h = hn::AverageRound(r3, r4);
        const auto use_src_v = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg_v, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, r1), threshold_v), hn::Ge(hn::AbsDiff(src, r2), threshold_v));
        const auto use_src_h = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg_h, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, r3), threshold_v), hn::Ge(hn::AbsDiff(src, r4), threshold_v));
        const auto new_v = hn::IfThenElse(use_src_v, src, avg_v);
        const auto new_h = hn::IfThenElse(use_src_h, src, avg_h);
        return hn::AverageRound(new_v, new_h);
    }
}

template <class DF, class VF>
HWY_INLINE VF process_mode6_vector(
    DF df,
    VF src_f,
    VF p1_f,
    VF p2_f,
    VF p3_f,
    VF p4_f,
    VF thresh_avg,
    VF thresh_max,
    VF thresh_mid
) {
    const auto quarter = hn::Set(df, 0.25f);
    const auto three = hn::Set(df, 3.0f);
    const auto eps = hn::Set(df, 1e-5f);
    const auto zero = hn::Zero(df);
    const auto one = hn::Set(df, 1.0f);

    const auto avg_refs_f = hn::Mul(hn::Add(hn::Add(p1_f, p2_f), hn::Add(p3_f, p4_f)), quarter);
    const auto diff_avg_src = hn::Sub(avg_refs_f, src_f);
    const auto avg_dif_f = hn::Abs(diff_avg_src);

    const auto d1 = hn::Abs(hn::Sub(p1_f, src_f));
    const auto d2 = hn::Abs(hn::Sub(p2_f, src_f));
    const auto d3 = hn::Abs(hn::Sub(p3_f, src_f));
    const auto d4 = hn::Abs(hn::Sub(p4_f, src_f));
    const auto max_dif = hn::Max(hn::Max(d1, d2), hn::Max(d3, d4));

    const auto two_src = hn::Add(src_f, src_f);
    const auto mid_dif_v = hn::Abs(hn::Sub(hn::Add(p1_f, p2_f), two_src));
    const auto mid_dif_h = hn::Abs(hn::Sub(hn::Add(p3_f, p4_f), two_src));

    const auto comp_avg = hn::Clamp(hn::Mul(three, hn::Sub(one, hn::Div(avg_dif_f, hn::Max(thresh_avg, eps)))), zero, one);
    const auto comp_max = hn::Clamp(hn::Mul(three, hn::Sub(one, hn::Div(max_dif, hn::Max(thresh_max, eps)))), zero, one);
    const auto comp_mid_v = hn::Clamp(hn::Mul(three, hn::Sub(one, hn::Div(mid_dif_v, hn::Max(thresh_mid, eps)))), zero, one);
    const auto comp_mid_h = hn::Clamp(hn::Mul(three, hn::Sub(one, hn::Div(mid_dif_h, hn::Max(thresh_mid, eps)))), zero, one);

    const auto product = hn::Mul(hn::Mul(comp_avg, comp_max), hn::Mul(comp_mid_v, comp_mid_h));

    // Vector pow(product, 0.1f) via Exp(0.1f * Log(Max(product, 1e-12f)))
    const auto log_val = hn::Log(df, hn::Max(product, hn::Set(df, 1e-12f)));
    const auto exp_val = hn::Exp(df, hn::Mul(hn::Set(df, 0.1f), log_val));
    const auto factor = hn::IfThenElse(hn::Gt(product, hn::Zero(df)), exp_val, hn::Zero(df));

    return hn::Add(src_f, hn::Mul(diff_avg_src, factor));
}

template <class DF, class VF>
HWY_INLINE VF process_mode7_vector(
    DF df,
    VF src_f,
    VF p1_f,
    VF p2_f,
    VF p3_f,
    VF p4_f,
    VF angle_org,
    VF angle_ref_h1,
    VF angle_ref_h2,
    VF angle_ref_w1,
    VF angle_ref_w2,
    VF thresh_avg,
    VF thresh_max,
    VF thresh_mid,
    float angle_boost,
    float max_angle
) {
    const auto max_diff1 = hn::Max(hn::Abs(hn::Sub(angle_ref_h1, angle_org)), hn::Abs(hn::Sub(angle_ref_h2, angle_org)));
    const auto max_diff2 = hn::Max(hn::Abs(hn::Sub(angle_ref_w1, angle_org)), hn::Abs(hn::Sub(angle_ref_w2, angle_org)));
    const auto max_angle_diff = hn::Max(max_diff1, max_diff2);

    const auto use_boost = hn::Le(max_angle_diff, hn::Set(df, max_angle));
    const auto boost_v = hn::Set(df, angle_boost);
    const auto t_avg = hn::IfThenElse(use_boost, hn::Mul(thresh_avg, boost_v), thresh_avg);
    const auto t_max = hn::IfThenElse(use_boost, hn::Mul(thresh_max, boost_v), thresh_max);
    const auto t_mid = hn::IfThenElse(use_boost, hn::Mul(thresh_mid, boost_v), thresh_mid);

    return process_mode6_vector(df, src_f, p1_f, p2_f, p3_f, p4_f, t_avg, t_max, t_mid);
}

template <int kSampleMode, bool kBlurFirst, class D32, class V32>
V32 process_reference_samples(
    D32 d32,
    V32 src,
    V32 r1,
    V32 r2,
    V32 r3,
    V32 r4,
    int threshold,
    int threshold1,
    int threshold2
) {
    const auto one = hn::Set(d32, 1);
    const auto threshold_v = hn::Set(d32, static_cast<std::int32_t>(threshold));

    if constexpr (kSampleMode == 1 || kSampleMode == 3) {
        const auto avg = hn::ShiftRight<1>(hn::Add(hn::Add(r1, r2), one));
        const auto use_src = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, r1), threshold_v), hn::Ge(hn::AbsDiff(src, r2), threshold_v));
        return hn::IfThenElse(use_src, src, avg);
    } else if constexpr (kSampleMode == 2) {
        auto avg1 = hn::ShiftRight<1>(hn::Add(hn::Add(r1, r2), one));
        const auto avg2 = hn::ShiftRight<1>(hn::Add(hn::Add(r3, r4), one));
        avg1 = hn::IfThenElse(hn::Gt(avg1, hn::Zero(d32)), hn::Sub(avg1, one), avg1);
        const auto avg = hn::ShiftRight<1>(hn::Add(hn::Add(avg1, avg2), one));
        const auto use_src = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg, src), threshold_v)
            : hn::Or(
                hn::Or(hn::Ge(hn::AbsDiff(r1, src), threshold_v), hn::Ge(hn::AbsDiff(r2, src), threshold_v)),
                hn::Or(hn::Ge(hn::AbsDiff(r3, src), threshold_v), hn::Ge(hn::AbsDiff(r4, src), threshold_v))
            );
        return hn::IfThenElse(use_src, src, avg);
    } else if constexpr (kSampleMode == 4) {
        const auto avg_v = hn::ShiftRight<1>(hn::Add(hn::Add(r1, r2), one));
        const auto avg_h = hn::ShiftRight<1>(hn::Add(hn::Add(r3, r4), one));
        const auto use_src_v = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg_v, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, r1), threshold_v), hn::Ge(hn::AbsDiff(src, r2), threshold_v));
        const auto use_src_h = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg_h, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, r3), threshold_v), hn::Ge(hn::AbsDiff(src, r4), threshold_v));
        const auto new_v = hn::IfThenElse(use_src_v, src, avg_v);
        const auto new_h = hn::IfThenElse(use_src_h, src, avg_h);
        return hn::ShiftRight<1>(hn::Add(hn::Add(new_v, new_h), one));
    } else if constexpr (kSampleMode == 5) {
        auto avg1 = hn::ShiftRight<1>(hn::Add(hn::Add(r1, r2), one));
        const auto avg2 = hn::ShiftRight<1>(hn::Add(hn::Add(r3, r4), one));
        avg1 = hn::IfThenElse(hn::Gt(avg1, hn::Zero(d32)), hn::Sub(avg1, one), avg1);
        const auto avg = hn::ShiftRight<1>(hn::Add(hn::Add(avg1, avg2), one));
        const auto threshold1_v = hn::Set(d32, static_cast<std::int32_t>(threshold1));
        const auto threshold2_v = hn::Set(d32, static_cast<std::int32_t>(threshold2));
        const auto max_diff = hn::Max(
            hn::AbsDiff(r1, src),
            hn::Max(hn::AbsDiff(r2, src), hn::Max(hn::AbsDiff(r3, src), hn::AbsDiff(r4, src)))
        );
        const auto use_src = hn::Or(
            hn::Or(hn::Ge(hn::AbsDiff(avg, src), threshold_v), hn::Ge(max_diff, threshold1_v)),
            hn::Or(
                hn::Ge(hn::AbsDiff(hn::Add(r1, r2), hn::ShiftLeft<1>(src)), threshold2_v),
                hn::Ge(hn::AbsDiff(hn::Add(r3, r4), hn::ShiftLeft<1>(src)), threshold2_v)
            )
        );
        return hn::IfThenElse(use_src, src, avg);
    } else if constexpr (kSampleMode == 6) {
        const hn::Rebind<float, D32> df;
        const auto src_f = hn::ConvertTo(df, src);
        const auto r1_f = hn::ConvertTo(df, r1);
        const auto r2_f = hn::ConvertTo(df, r2);
        const auto r3_f = hn::ConvertTo(df, r3);
        const auto r4_f = hn::ConvertTo(df, r4);
        const auto thresh_avg = hn::Set(df, static_cast<float>(threshold));
        const auto thresh_max = hn::Set(df, static_cast<float>(threshold1));
        const auto thresh_mid = hn::Set(df, static_cast<float>(threshold2));

        const auto blended_f = process_mode6_vector(df, src_f, r1_f, r2_f, r3_f, r4_f, thresh_avg, thresh_max, thresh_mid);
        return hn::NearestInt(blended_f);
    } else {
        static_assert(kSampleMode >= 1 && kSampleMode <= 7);
        return src;
    }
}
