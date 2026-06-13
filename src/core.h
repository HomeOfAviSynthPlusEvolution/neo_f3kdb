#pragma once

#include "aligned_buffer.hpp"
#include <cstdint>
#include <vector>

using neo_f3kdb::AlignedBuffer;

// Lightweight replacement for old dualsynth1 headers to avoid conflict with dualsynth2
struct DSFormat {
  bool IsFamilyYUV {true};
  int SSW {0};
  int SSH {0};
  int BitsPerSample {8};
  int Planes {3};
};

struct DSVideoInfo {
  DSFormat Format;
  int Width {0};
  int Height {0};
  int Frames {0};
};

#include "f3kdb.h"
#include "process_plane_context.h"
#include "compiler_compat.h"

typedef struct _pixel_dither_info {
    alignas(4) std::int8_t ref1;
    std::int8_t ref2;
    std::int16_t change;
} pixel_dither_info;

static_assert(sizeof(pixel_dither_info) == 4, "Something wrong in pixel_dither_info");

typedef struct _process_plane_params
{
    const unsigned char *src_plane_ptr;
    int src_pitch;

    unsigned char *dst_plane_ptr;
    int dst_pitch;

    int plane_width_in_pixels;
    int plane_height_in_pixels;

    PIXEL_MODE input_mode;
    int input_depth;
    PIXEL_MODE output_mode;
    int output_depth;

    std::uint16_t threshold;
    std::uint16_t threshold1;
    std::uint16_t threshold2;
    float angle_boost;
    float max_angle;
    pixel_dither_info *info_ptr_base;
    int info_stride;
    
    std::int16_t* grain_buffer;
    int grain_buffer_stride;

    int plane;

    unsigned char width_subsampling;
    unsigned char height_subsampling;
    
    int pixel_max;
    int pixel_min;
    
    // Helper functions
    inline int get_dst_width() const {
        return output_mode == HIGH_BIT_DEPTH_INTERLEAVED ? plane_width_in_pixels * 2 : plane_width_in_pixels;
    }
    inline int get_dst_height() const {
        return plane_height_in_pixels;
    }
    inline int get_src_width() const {
        return input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? plane_width_in_pixels * 2 : plane_width_in_pixels;
    }
    inline int get_src_height() const {
        return plane_height_in_pixels;
    }
} process_plane_params;

typedef void (*process_plane_impl_t)(const process_plane_params& params, process_plane_context* context);

class f3kdb_core_t {
private:
    process_plane_impl_t _process_plane_impl;
        
    AlignedBuffer<pixel_dither_info, 128> _y_info;
    AlignedBuffer<pixel_dither_info, 128> _cb_info;
    AlignedBuffer<pixel_dither_info, 128> _cr_info;
    
    process_plane_context _y_context;
    process_plane_context _cb_context;
    process_plane_context _cr_context;
    
    AlignedBuffer<std::int16_t, 128> _grain_buffer_y;
    AlignedBuffer<std::int16_t, 128> _grain_buffer_c;

    std::vector<std::int32_t> _grain_buffer_offsets;

    DSVideoInfo _video_info;
    f3kdb_params_t _params;

    OPTIMIZATION_MODE _opt;

    void init(void);
    void init_frame_luts(void);

    void destroy_frame_luts(void);

    f3kdb_core_t(const f3kdb_core_t&);
    f3kdb_core_t operator=(const f3kdb_core_t&);
    
public:
    f3kdb_core_t(DSVideoInfo vi, const f3kdb_params_t params, OPTIMIZATION_MODE opt);
    virtual ~f3kdb_core_t();

    void process_plane(int frame_index, int plane, unsigned char* dst_frame_ptr, int dst_pitch, const unsigned char* src_frame_ptr, int src_pitch);
};
