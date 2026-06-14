#ifndef PIXEL_PROC_C_HIGH_F_S_DITHERING_H
#define PIXEL_PROC_C_HIGH_F_S_DITHERING_H

#include <cstdlib>
#include <cstring>
#include <new>

#include "core/constants.hpp"
#include "core/kernel.hpp"
#include "core/pixel_proc_common.hpp"

namespace neo_f3kdb::core::dither {
    
// #define DUMP_DATA

class FloydSteinbergDither {
public:
    using ErrorType = unsigned short;

    FloydSteinbergDither(char* inline_buffer, int inline_buffer_size, int frame_width, int output_depth)
        : output_depth_(output_depth),
          row_pitch_(frame_width + 2),
          frame_width_(frame_width)
    {
        const int size_needed = (frame_width + 2) * 2 * static_cast<int>(sizeof(ErrorType));
        if (inline_buffer_size < size_needed)
        {
            error_buffer_ = static_cast<ErrorType*>(std::malloc(size_needed));
            buffer_needs_dealloc_ = true;
        } else {
            error_buffer_ = reinterpret_cast<ErrorType*>(inline_buffer);
        }
        std::memset(error_buffer_, 0, size_needed);
        current_px_error_ = error_buffer_ + 1;

#ifdef DUMP_DATA
        char file_name[256];
        sprintf(file_name, "fsdither_dump_stage0_%d", frame_width);
        debug_dump_fd_[0] = fopen(file_name, "wb");
        sprintf(file_name, "fsdither_dump_stage1_%d", frame_width);
        debug_dump_fd_[1] = fopen(file_name, "wb");
        sprintf(file_name, "fsdither_dump_stage2_%d", frame_width);
        debug_dump_fd_[2] = fopen(file_name, "wb");
#endif
    }

    ~FloydSteinbergDither()
    {
        if (buffer_needs_dealloc_)
        {
            std::free(error_buffer_);
            error_buffer_ = nullptr;
        }
#ifdef DUMP_DATA
        for (int i = 0; i < sizeof(debug_dump_fd_) / sizeof(FILE*); i++)
        {
            if (debug_dump_fd_[i])
            {
                fclose(debug_dump_fd_[i]);
            }
        }
#endif
    }

    void next_pixel()
    {
        current_px_error_++;
        processed_pixels_in_current_line_++;
    }

    void next_row()
    {
        row_pitch_ = -row_pitch_;
        current_px_error_ = error_buffer_ + (row_pitch_ >> 31) * row_pitch_;
        std::memset(current_px_error_ + row_pitch_, 0, std::abs(row_pitch_) * sizeof(ErrorType));
        current_px_error_++;
        processed_pixels_in_current_line_ = 0;
    }

    int dither(int pixel)
    {
        if (processed_pixels_in_current_line_ >= frame_width_)
        {
            // outside plane, can occur in SSE code
            return pixel;
        }
#ifndef FS_DITHER_SKIP_PRE_CLAMP
        pixel = clamp_pixel(pixel, kPixelMin, kPixelMax);
#endif
#ifdef DUMP_DATA
        fwrite(&pixel, 4, 1, debug_dump_fd_[0]);
#endif
        pixel += *current_px_error_;
#ifdef DUMP_DATA
        fwrite(&pixel, 4, 1, debug_dump_fd_[1]);
#endif
        pixel = clamp_pixel(pixel, kPixelMin, kPixelMax);
#ifdef DUMP_DATA
        fwrite(&pixel, 4, 1, debug_dump_fd_[2]);
#endif
        const int new_error = pixel & ((1 << (INTERNAL_BIT_DEPTH - output_depth_)) - 1);
        *(current_px_error_ + 1) += (new_error * 7) >> 4;
        *(current_px_error_ + row_pitch_ - 1) += (new_error * 3) >> 4;
        *(current_px_error_ + row_pitch_) += (new_error * 5) >> 4;
        *(current_px_error_ + row_pitch_ + 1) += (new_error * 1) >> 4;
        return pixel;
    }

private:
    static constexpr int kPixelMax = (1 << INTERNAL_BIT_DEPTH) - 1;
    static constexpr int kPixelMin = 0;

    int output_depth_ = 0;
    ErrorType* error_buffer_ = nullptr;
    bool buffer_needs_dealloc_ = false;
    ErrorType* current_px_error_ = nullptr;
    int row_pitch_ = 0;
    int frame_width_ = 0;
    int processed_pixels_in_current_line_ = 0;
#ifdef DUMP_DATA
    FILE* debug_dump_fd_[3] {};
#endif
};

} // namespace neo_f3kdb::core::dither

namespace pixel_proc_high_f_s_dithering {

    using FloydSteinbergDither = neo_f3kdb::core::dither::FloydSteinbergDither;

    static inline FloydSteinbergDither& state(void* context)
    {
        return *reinterpret_cast<FloydSteinbergDither*>(context);
    }

    static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth)
    {
        new (context_buffer) FloydSteinbergDither(
            context_buffer + sizeof(FloydSteinbergDither),
            CONTEXT_BUFFER_SIZE - static_cast<int>(sizeof(FloydSteinbergDither)),
            frame_width,
            output_depth
        );
    }

    static inline void destroy_context(void* context)
    {
        state(context).~FloydSteinbergDither();
    }

    static __forceinline void next_pixel(void* context)
    {
        state(context).next_pixel();
    }

    static __forceinline void next_row(void* context)
    {
        state(context).next_row();
    }

    static __forceinline int dither(void* context, int pixel, int row, int column)
    {
        return state(context).dither(pixel);
    }

    static inline int upsample(void* context, unsigned char pixel)
    {
        return neo_f3kdb::core::pixel_proc_common::upsample(context, pixel);
    }

    static inline int downsample(void* context, int pixel, int row, int column, int pixel_min, int pixel_max, int output_depth)
    {
        pixel = dither(context, pixel, row, column);
        return neo_f3kdb::core::pixel_proc_common::downsample(pixel, pixel_min, pixel_max, output_depth);
    }

    static inline int avg_2(void* context, int pixel1, int pixel2)
    {
        return neo_f3kdb::core::pixel_proc_common::avg_2(context, pixel1, pixel2);
    }

    static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4)
    {
        return neo_f3kdb::core::pixel_proc_common::avg_4(context, pixel1, pixel2, pixel3, pixel4);
    }

};

#endif // PIXEL_PROC_C_HIGH_F_S_DITHERING_H
