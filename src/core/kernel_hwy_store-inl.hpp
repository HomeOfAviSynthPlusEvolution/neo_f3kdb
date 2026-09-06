// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <class D32, class Pixel, class V32>
void store_pixels_from_u32(D32 d32, Pixel* dest_ptr, V32 values, bool is_chroma = false) {
    if constexpr (std::is_same_v<Pixel, std::uint8_t>) {
        const hn::Rebind<std::uint16_t, D32> d16;
        const hn::Rebind<std::uint8_t, D32> d8;
        const auto out16 = hn::DemoteTo(d16, values);
        const auto out8 = hn::DemoteTo(d8, out16);
        hn::StoreU(out8, d8, dest_ptr);
    } else if constexpr (std::is_same_v<Pixel, float>) {
        const hn::Rebind<float, D32> df;
        auto out_f = hn::Mul(hn::ConvertTo(df, values), hn::Set(df, 1.0f / 65535.0f));
        if (is_chroma) {
            out_f = hn::Sub(out_f, hn::Set(df, 0.5f));
        }
        hn::StoreU(out_f, df, dest_ptr);
    } else {
        const hn::Rebind<std::uint16_t, D32> d16;
        const auto out16 = hn::DemoteTo(d16, values);
        hn::StoreU(out16, d16, dest_ptr);
    }
}
