// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <int kSampleMode, bool kBlurFirst, class D32, class V32>
V32 process_reference_samples(
    D32 d32,
    V32 src,
    const std::int32_t* ref1,
    const std::int32_t* ref2,
    const std::int32_t* ref3,
    const std::int32_t* ref4,
    int threshold,
    int threshold1,
    int threshold2
) {
    const auto one = hn::Set(d32, 1);
    const auto threshold_v = hn::Set(d32, static_cast<std::int32_t>(threshold));

    if constexpr (kSampleMode == 1 || kSampleMode == 3) {
        const auto r1 = hn::LoadU(d32, ref1);
        const auto r2 = hn::LoadU(d32, ref2);
        const auto avg = hn::ShiftRight<1>(hn::Add(hn::Add(r1, r2), one));
        const auto use_src = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, r1), threshold_v), hn::Ge(hn::AbsDiff(src, r2), threshold_v));
        return hn::IfThenElse(use_src, src, avg);
    } else if constexpr (kSampleMode == 2) {
        const auto r1 = hn::LoadU(d32, ref1);
        const auto r2 = hn::LoadU(d32, ref2);
        const auto r3 = hn::LoadU(d32, ref3);
        const auto r4 = hn::LoadU(d32, ref4);
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
        const auto ref_v1 = hn::LoadU(d32, ref1);
        const auto ref_v2 = hn::LoadU(d32, ref2);
        const auto ref_h1 = hn::LoadU(d32, ref3);
        const auto ref_h2 = hn::LoadU(d32, ref4);
        const auto avg_v = hn::ShiftRight<1>(hn::Add(hn::Add(ref_v1, ref_v2), one));
        const auto avg_h = hn::ShiftRight<1>(hn::Add(hn::Add(ref_h1, ref_h2), one));
        const auto use_src_v = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg_v, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, ref_v1), threshold_v), hn::Ge(hn::AbsDiff(src, ref_v2), threshold_v));
        const auto use_src_h = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg_h, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, ref_h1), threshold_v), hn::Ge(hn::AbsDiff(src, ref_h2), threshold_v));
        const auto new_v = hn::IfThenElse(use_src_v, src, avg_v);
        const auto new_h = hn::IfThenElse(use_src_h, src, avg_h);
        return hn::ShiftRight<1>(hn::Add(hn::Add(new_v, new_h), one));
    } else {
        static_assert(kSampleMode == 5);
        const auto r1 = hn::LoadU(d32, ref1);
        const auto r2 = hn::LoadU(d32, ref2);
        const auto r3 = hn::LoadU(d32, ref3);
        const auto r4 = hn::LoadU(d32, ref4);
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
