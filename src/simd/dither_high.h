
#include "../impl_dispatch.h"
#include "../compiler_compat.h"

#define FS_DITHER_SKIP_PRE_CLAMP

#include "../pixel_proc_c_high_f_s_dithering.h"
#include "../pixel_proc_c_high_ordered_dithering.h"

#include <assert.h>
#include <emmintrin.h>

namespace dither_high
{
    static SIMD::data_type _ordered_dithering_threshold_map[16] [2];
    static SIMD::data_type _ordered_dithering_threshold_map_yuy2[16] [8];
    static volatile bool _threshold_map_initialized = false;

    static __inline void init_ordered_dithering()
    {
        if (!_threshold_map_initialized) {
            SIMD::data_type threhold_row;
            SIMD::data_type zero = SIMD::_setzero();
            for (int i = 0; i < 16; i++) 
            {
                threhold_row = *(SIMD::data_type*)pixel_proc_high_ordered_dithering::THRESHOLD_MAP[i];
                    
                SIMD::data_type part_0 = SIMD::_unpacklo_epi8(threhold_row, zero);
                SIMD::data_type part_1 = SIMD::_unpackhi_epi8(threhold_row, zero);

                if (INTERNAL_BIT_DEPTH < 16)
                {
                    part_0 = SIMD::_srli_epi16(part_0, 16 - INTERNAL_BIT_DEPTH);
                    part_1 = SIMD::_srli_epi16(part_1, 16 - INTERNAL_BIT_DEPTH);
                }
                _ordered_dithering_threshold_map[i][0] = part_0;
                _ordered_dithering_threshold_map[i][1] = part_1;
                
                SIMD::data_type tmp = SIMD::_unpacklo_epi8(part_0, part_0);
                _ordered_dithering_threshold_map_yuy2[i][0] = SIMD::_unpacklo_epi16(part_0, tmp);
                _ordered_dithering_threshold_map_yuy2[i][1] = SIMD::_unpackhi_epi16(part_0, tmp);

                tmp = SIMD::_unpackhi_epi8(part_0, part_0);
                _ordered_dithering_threshold_map_yuy2[i][2] = SIMD::_unpacklo_epi16(part_1, tmp);
                _ordered_dithering_threshold_map_yuy2[i][3] = SIMD::_unpackhi_epi16(part_1, tmp);

                tmp = SIMD::_unpacklo_epi8(part_1, part_1);
                _ordered_dithering_threshold_map_yuy2[i][4] = SIMD::_unpacklo_epi16(part_0, tmp);
                _ordered_dithering_threshold_map_yuy2[i][5] = SIMD::_unpackhi_epi16(part_0, tmp);

                tmp = SIMD::_unpackhi_epi8(part_1, part_1);
                _ordered_dithering_threshold_map_yuy2[i][6] = SIMD::_unpacklo_epi16(part_1, tmp);
                _ordered_dithering_threshold_map_yuy2[i][7] = SIMD::_unpackhi_epi16(part_1, tmp);
            }
            _mm_mfence();
            _threshold_map_initialized = true;
        }
    }

    static void init_ordered_dithering_with_output_depth(char context_buffer[CONTEXT_BUFFER_SIZE], int output_depth)
    {
        assert(_threshold_map_initialized);

        SIMD::count_type shift = SIMD::count_set(output_depth - 8);

        for (int i = 0; i < 16; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                SIMD::data_type item = _ordered_dithering_threshold_map[i][j];
                item = SIMD::_srl_epi16(item, shift);
                SIMD::_store((SIMD::data_type*)(context_buffer + (i * 2 + j) * 16), item);
            }
        }
    }

    template <int dither_algo>
    static __inline void init(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) 
    {
        if (dither_algo == DA_HIGH_FLOYD_STEINBERG_DITHERING)
        {
            pixel_proc_high_f_s_dithering::init_context(context_buffer, frame_width, output_depth);
        } else if (dither_algo == DA_HIGH_ORDERED_DITHERING) {
            init_ordered_dithering();
            init_ordered_dithering_with_output_depth(context_buffer, output_depth);
        }
    }

    template <int dither_algo>
    static __inline void complete(void* context) 
    {
        if (dither_algo == DA_HIGH_FLOYD_STEINBERG_DITHERING)
        {
            pixel_proc_high_f_s_dithering::destroy_context(context);
        }
    }
    
    template <int dither_algo>
    static __forceinline SIMD::data_type dither(void* context, SIMD::data_type pixels, int row, int column)
    {
        switch (dither_algo)
        {
        case DA_HIGH_NO_DITHERING:
            return pixels;
        case DA_HIGH_ORDERED_DITHERING:
            {
            // row: use lowest 4 bits as index, mask = 0b00001111 = 15
            // column: always multiples of 8, so use 8 (bit 4) as selector, mask = 0b00001000
            assert((column & 7) == 0);
            SIMD::data_type threshold = SIMD::_load((SIMD::data_type*)((char*)context + ( ( (row & 15) * 2 ) + ( (column & 8) >> 3 ) ) * 16 ) );
            return SIMD::_adds_epu16(pixels, threshold);
            }
        case DA_HIGH_FLOYD_STEINBERG_DITHERING:
            // fixme, remove shitty compat
            // due to an ICC bug, accessing pixels using union will give us incorrect results
            // so we have to use a buffer here
            // tested on ICC 12.0.1024.2010
            alignas(SIMD::align) unsigned short buffer[SIMD::width_16];
            SIMD::_store((SIMD::data_type*)buffer, pixels);
            for (int i = 0; i < SIMD::width_16; i++)
            {
                buffer[i] = (unsigned short)pixel_proc_high_f_s_dithering::dither(context, buffer[i], row, column + i);
                pixel_proc_high_f_s_dithering::next_pixel(context);
            }
            return SIMD::_load((SIMD::data_type*)buffer);
        case DA_16BIT_INTERLEAVED:
            return SIMD::_setzero();
            break;
        default:
            abort();
            return SIMD::_setzero();
        }
    }
    
    template <int dither_algo>
    static __inline void next_row(void* context)
    {
        if (dither_algo == DA_HIGH_FLOYD_STEINBERG_DITHERING)
        {
            pixel_proc_high_f_s_dithering::next_row(context);
        }
    }
};
