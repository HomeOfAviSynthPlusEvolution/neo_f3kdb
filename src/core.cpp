#include <stdlib.h>
#include <stdarg.h>
#include <memory.h>
#include <assert.h>

#include "core.h"
#include "constants.h"
#include "random.h"
#include "impl_dispatch.h"

#ifdef _WIN32
#include <intrin.h>
#endif

void f3kdb_core_t::destroy_frame_luts(void)
{
    _aligned_free(_y_info);
    _aligned_free(_cb_info);
    _aligned_free(_cr_info);
    
    _y_info = NULL;
    _cb_info = NULL;
    _cr_info = NULL;
    
    _aligned_free(_grain_buffer_y);
    _aligned_free(_grain_buffer_c);
    
    _grain_buffer_y = NULL;
    _grain_buffer_c = NULL;

    free(_grain_buffer_offsets);
    _grain_buffer_offsets = NULL;
    
    // contexts are likely to be dependent on lut, so they must also be destroyed
    destroy_context(&_y_context);
    destroy_context(&_cb_context);
    destroy_context(&_cr_context);
}

static int inline min_multi( int first, ... )
{
    int ret = first, i = first;
    va_list marker;

    va_start( marker, first );
    while( i >= 0 )
    {
        if (i < ret)
        {
            ret = i;
        }
        i = va_arg( marker, int);
    }
    va_end( marker );
    return ret;
}

static int get_frame_lut_stride(int width_in_pixels)
{
    // whole multiples of alignment, so SSE codes don't need to check boundaries
    int width = width_in_pixels;
    return (((width - 1) | (FRAME_LUT_ALIGNMENT - 1)) + 1);
}

static short* generate_grain_buffer(size_t item_count, RANDOM_ALGORITHM algo, int& seed, double param, int range)
{
    short* buffer = (short*)_aligned_malloc(item_count * sizeof(short), FRAME_LUT_ALIGNMENT);
    for (size_t i = 0; i < item_count; i++)
    {
        *(buffer + i) = random(algo, seed, range, param);
    }
    return buffer;
}

void f3kdb_core_t::init_frame_luts(void)
{
    destroy_frame_luts();

    int seed = 0x92D68CA2 - _params.seed;

    seed ^= (_video_info.Width << 16) ^ _video_info.Height;
    seed ^= (_video_info.Frames << 16) ^ _video_info.Frames;

    int height_in_pixels = _video_info.Height;
    int width_in_pixels =  _video_info.Width;

    int y_stride;
    y_stride = get_frame_lut_stride(width_in_pixels);

    int y_size = sizeof(pixel_dither_info) * y_stride * height_in_pixels;
    _y_info = (pixel_dither_info*)_aligned_malloc(y_size, FRAME_LUT_ALIGNMENT);

    // ensure unused items are also initialized
    memset(_y_info, 0, y_size);

    int c_stride;
    c_stride = get_frame_lut_stride(width_in_pixels >> _video_info.Format.SSW);
    int c_size = sizeof(pixel_dither_info) * c_stride * (height_in_pixels >> _video_info.Format.SSH);
    _cb_info = (pixel_dither_info*)_aligned_malloc(c_size, FRAME_LUT_ALIGNMENT);
    _cr_info = (pixel_dither_info*)_aligned_malloc(c_size, FRAME_LUT_ALIGNMENT);

    memset(_cb_info, 0, c_size);
    memset(_cr_info, 0, c_size);

    pixel_dither_info *y_info_ptr, *cb_info_ptr, *cr_info_ptr;

    int width_subsamp = _video_info.Format.SSW;
    int height_subsamp = _video_info.Format.SSH;

    for (int y = 0; y < height_in_pixels; y++)
    {
        y_info_ptr = _y_info + y * y_stride;
        cb_info_ptr = _cb_info + (y >> height_subsamp) * c_stride;
        cr_info_ptr = _cr_info + (y >> height_subsamp) * c_stride;

        for (int x = 0; x < width_in_pixels; x++)
        {
            pixel_dither_info info_y = {0, 0, 0};
            info_y.change = random(_params.random_algo_grain, seed, _params.grainY, _params.random_param_grain);

            int x_range = min_multi(_params.range, x, width_in_pixels - x - 1, -1);
            int y_range = min_multi(_params.range, y, height_in_pixels - y - 1, -1);
            int cur_range;
            switch (_params.sample_mode)
            {
            case 1:
                cur_range = y_range;
                break;
            
            case 3:
                cur_range = x_range;
                break;
            
            case 2:
            case 4:
                cur_range = min_multi(x_range, y_range, -1);
                break;

            default:
                break;
            }

            if (cur_range > 0) {
                info_y.ref1 = (signed char)random(_params.random_algo_ref, seed, cur_range, _params.random_param_ref);
                if (_params.sample_mode == 2)
                {
                    info_y.ref2 = (signed char)random(_params.random_algo_ref, seed, cur_range, _params.random_param_ref);
                }
                if (_params.sample_mode > 0)
                {
                    info_y.ref1 = abs(info_y.ref1);
                    info_y.ref2 = abs(info_y.ref2);
                }
            }

            *y_info_ptr = info_y;

            bool should_set_c = false;
            should_set_c = ((x & ( ( 1 << width_subsamp ) - 1)) == 0 && 
                (y & ( ( 1 << height_subsamp ) - 1)) == 0);

            if (should_set_c) {
                pixel_dither_info info_cb = info_y;
                pixel_dither_info info_cr = info_cb;

                // don't shift ref values here, since subsampling of width and height may be different
                // shift them in actual processing

                info_cb.change = random(_params.random_algo_grain, seed, _params.grainC, _params.random_param_grain);
                info_cr.change = random(_params.random_algo_grain, seed, _params.grainC, _params.random_param_grain);

                *cb_info_ptr = info_cb;
                *cr_info_ptr = info_cr;
                cb_info_ptr++;
                cr_info_ptr++;
            }
            y_info_ptr++;
        }
    }

    int multiplier = _params.dynamic_grain ? 3 : 1;
    int item_count = width_in_pixels;

    // add some safety margin and align it
    item_count += 255;
    item_count &= 0xffffff80;

    item_count *= height_in_pixels;

    _grain_buffer_y = generate_grain_buffer(
        item_count * multiplier,
        _params.random_algo_grain,
        seed,
        _params.random_param_grain,
        _params.grainY);

    // we always generate a full-sized buffer to simplify offset calculation
    _grain_buffer_c = generate_grain_buffer(
        item_count * multiplier,
        _params.random_algo_grain,
        seed,
        _params.random_param_grain,
        _params.grainC);

    if (_params.dynamic_grain)
    {
        // Pre-generate offset here so that result is deterministic even if we request frame in different order
        _grain_buffer_offsets = (int*)malloc(sizeof(int) * _video_info.Frames);
        for (int i = 0; i < _video_info.Frames; i++)
        {
            int offset = item_count + random(RANDOM_ALGORITHM_UNIFORM, seed, item_count, DEFAULT_RANDOM_PARAM);
            offset &= 0xfffffff0; // align to 16-byte for SSE codes

            assert(offset >= 0);

            _grain_buffer_offsets[i] = offset;
        }
    }
}

f3kdb_core_t::f3kdb_core_t(DSVideoInfo vi, const f3kdb_params_t params, OPTIMIZATION_MODE opt) :
    _process_plane_impl(NULL),
    _y_info(NULL),
    _cb_info(NULL),
    _cr_info(NULL),
    _grain_buffer_y(NULL),
    _grain_buffer_c(NULL),
    _grain_buffer_offsets(NULL),
    _video_info(vi),
    _opt(opt),
    _params(params)
{
    this->init();
}

f3kdb_core_t::~f3kdb_core_t()
{
    destroy_frame_luts();
}

static __inline int select_impl_index(int sample_mode, bool blur_first)
{
    assert(sample_mode != 0);
    return sample_mode * 2 + (blur_first ? 0 : 1) - 1;
}

void f3kdb_core_t::init(void) 
{
    init_context(&_y_context);
    init_context(&_cb_context);
    init_context(&_cr_context);

    init_frame_luts();

    const process_plane_impl_t* impl_table = process_plane_impls[_params.dither_algo][(int)_opt];
    _process_plane_impl = impl_table[select_impl_index(_params.sample_mode, _params.blur_first)];
}

void f3kdb_core_t::process_plane(int frame_index, int plane, unsigned char* dst_frame_ptr, int dst_pitch, const unsigned char* src_frame_ptr, int src_pitch)
{
    process_plane_params params;

    memset(&params, 0, sizeof(process_plane_params));

    params.src_plane_ptr = src_frame_ptr;
    params.src_pitch = src_pitch;

    params.dst_plane_ptr = dst_frame_ptr;
    params.dst_pitch = dst_pitch;

    params.input_mode = _video_info.Format.BitsPerSample == 8 ? LOW_BIT_DEPTH : HIGH_BIT_DEPTH_INTERLEAVED;
    params.input_depth = _video_info.Format.BitsPerSample;
    params.output_mode = _params.output_depth <= 8 ? LOW_BIT_DEPTH : HIGH_BIT_DEPTH_INTERLEAVED;
    params.output_depth = _params.output_depth;

    params.plane = plane;
    
    params.width_subsampling = plane == 0 ? 0 : _video_info.Format.SSW;
    params.height_subsampling = plane == 0 ? 0 : _video_info.Format.SSH;

    params.plane_width_in_pixels = plane == 0 ? _video_info.Width : (_video_info.Width >> _video_info.Format.SSW);
    params.plane_height_in_pixels = plane == 0 ? _video_info.Height : (_video_info.Height >> _video_info.Format.SSH);

    params.info_stride = get_frame_lut_stride(params.plane_width_in_pixels);
    params.grain_buffer_stride = get_frame_lut_stride(params.plane_width_in_pixels);

    process_plane_context* context;

    int grain_setting = 0;

    switch (plane)
    {
    case 0:
        params.info_ptr_base = _y_info;
        params.threshold = _params.Y;
        params.pixel_max = _params.keep_tv_range ? TV_RANGE_Y_MAX : FULL_RANGE_Y_MAX;
        params.pixel_min = _params.keep_tv_range ? TV_RANGE_Y_MIN : FULL_RANGE_Y_MIN;
        params.grain_buffer = _grain_buffer_y;
        grain_setting = _params.grainY;
        context = &_y_context;
        break;
    case 1:
        params.info_ptr_base = _cb_info;
        params.threshold = _params.Cb;
        params.pixel_max = _params.keep_tv_range ? TV_RANGE_C_MAX : FULL_RANGE_C_MAX;
        params.pixel_min = _params.keep_tv_range ? TV_RANGE_C_MIN : FULL_RANGE_C_MIN;
        params.grain_buffer = _grain_buffer_c;
        grain_setting = _params.grainC;
        context = &_cb_context;
        break;
    case 2:
        params.info_ptr_base = _cr_info;
        params.threshold = _params.Cr;
        params.pixel_max = _params.keep_tv_range ? TV_RANGE_C_MAX : FULL_RANGE_C_MAX;
        params.pixel_min = _params.keep_tv_range ? TV_RANGE_C_MIN : FULL_RANGE_C_MIN;
        params.grain_buffer = _grain_buffer_c;
        grain_setting = _params.grainC;
        context = &_cr_context;
        break;
    default:
        abort();
    }
    
    if (_grain_buffer_offsets)
    {
        params.grain_buffer += _grain_buffer_offsets[frame_index % _video_info.Frames];
    }

    bool copy_plane = false;
    if (_video_info.Format.BitsPerSample == _params.output_depth &&
        grain_setting == 0 &&
        params.threshold == 0)
    {
        copy_plane = true;
    }

    if (copy_plane) {
        // no need to process
        int line_size = params.get_src_width();
        auto src = src_frame_ptr;
        auto dst = dst_frame_ptr;
        if (line_size == src_pitch && src_pitch == dst_pitch)
        {
            memcpy(dst, src, line_size * params.get_src_height());
        } else {
            for (int row = 0; row < params.get_src_height(); row++) 
            {
                memcpy(dst, src, line_size);
                src += src_pitch;
                dst += dst_pitch;
            }
        }
        return;
    }

    _process_plane_impl(params, context);
}
