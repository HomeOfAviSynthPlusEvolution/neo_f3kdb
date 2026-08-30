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
    } else {
        static_assert(kSampleMode == 5);
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
    }
}
