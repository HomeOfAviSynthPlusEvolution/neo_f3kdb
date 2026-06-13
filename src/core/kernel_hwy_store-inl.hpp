// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <class D32, class Pixel, class V32>
void store_pixels_from_u32(D32 d32, Pixel* dest_ptr, V32 values) {
    if constexpr (std::is_same_v<Pixel, std::uint8_t>) {
        const hn::Rebind<std::uint16_t, D32> d16;
        const hn::Rebind<std::uint8_t, D32> d8;
        const auto out16 = hn::DemoteTo(d16, values);
        const auto out8 = hn::DemoteTo(d8, out16);
        hn::StoreU(out8, d8, dest_ptr);
    } else {
        const hn::Rebind<std::uint16_t, D32> d16;
        const auto out16 = hn::DemoteTo(d16, values);
        hn::StoreU(out16, d16, dest_ptr);
    }
}
