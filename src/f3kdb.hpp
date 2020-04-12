/*
 * Copyright 2020 Xinyue Lu
 *
 * Temporal Median - filter.
 *
 */

#pragma once

#include <execution>

#include "compiler_compat.h"
#include "core.h"
#include "constants.h"
#include "impl_dispatch.h"

int GetCPUFlags();

struct F3KDB final : Filter {
  f3kdb_params_t ep;
  f3kdb_core_t* engine;
  InDelegator* _in;
  bool crop;
  char error_msg[1024];
  DSVideoInfo out_vi;

  const char* VSName() const override { return "Deband"; }
  const char* AVSName() const override { return "neo_f3kdb"; }
  const MtMode AVSMode() const override { return MT_NICE_FILTER; }
  const VSFilterMode VSMode() const override { return fmParallel; }
  const std::vector<Param> Params() const override {
    return std::vector<Param> {
      Param {"clip", Clip, false, true, true, false},
      Param {"range", Integer},
      Param {"y", Integer},
      Param {"cb", Integer},
      Param {"cr", Integer},
      Param {"grainy", Integer},
      Param {"grainc", Integer},
      Param {"sample_mode", Integer},
      Param {"seed", Integer},
      Param {"blur_first", Boolean},
      Param {"dynamic_grain", Boolean},
      Param {"dither_algo", Integer},
      Param {"keep_tv_range", Boolean},
      Param {"output_depth", Integer},
      Param {"random_algo_ref", Integer},
      Param {"random_algo_grain", Integer},
      Param {"random_param_ref", Float},
      Param {"random_param_grain", Float},
      Param {"preset", String}
    };
  }
  void Initialize(InDelegator* in, DSVideoInfo in_vi, FetchFrameFunctor* fetch_frame) override
  {
    Filter::Initialize(in, in_vi, fetch_frame);
    std::string preset;
    in->Read("preset", preset);
    if (preset == "depth")
      ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = 0;
    else if (preset == "low")
      ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = 32;
    else if (preset == "medium")
      ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = 48;
    else if (preset == "high")
      ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = 64;
    else if (preset == "veryhigh")
      ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = 80;
    else if (preset == "nograin")
      ep.grainY = ep.grainC = 0;
    else if (preset == "luma")
      ep.Cb = ep.Cr = ep.grainC = 0;
    else if (preset == "chroma")
      ep.Y = ep.grainY = 0;
    int tmp;
    in->Read("range", ep.range);
    in->Read("Y", ep.Y);
    in->Read("Cb", ep.Cb);
    in->Read("Cr", ep.Cr);
    in->Read("grainY", ep.grainY);
    in->Read("grainC", ep.grainC);
    in->Read("sample_mode", ep.sample_mode);
    in->Read("seed", ep.seed);
    in->Read("blur_first", ep.blur_first);
    in->Read("dynamic_grain", ep.dynamic_grain);
    tmp = static_cast<int>(ep.dither_algo);
    in->Read("dither_algo", tmp);
    ep.dither_algo = static_cast<DITHER_ALGORITHM>(tmp);
    in->Read("keep_tv_range", ep.keep_tv_range);
    in->Read("output_depth", ep.output_depth);
    tmp = static_cast<int>(ep.random_algo_ref);
    in->Read("random_algo_ref", tmp);
    ep.random_algo_ref = static_cast<RANDOM_ALGORITHM>(tmp);
    tmp = static_cast<int>(ep.random_algo_grain);
    in->Read("random_algo_grain", tmp);
    ep.random_algo_grain = static_cast<RANDOM_ALGORITHM>(tmp);
    in->Read("random_param_ref", ep.random_param_ref);
    in->Read("random_param_grain", ep.random_param_grain);

    OPTIMIZATION_MODE opt = IMPL_C;
    int CPUFlags = GetCPUFlags();

    if (CPUFlags & CPUF_SSE4_1)
      opt = IMPL_SSE4;

    #define INVALID_PARAM_IF(cond) \
    do { if (cond) { throw("Invalid parameter condition: " #cond); } } while (0)

    INVALID_PARAM_IF(in_vi.Format.IsFamilyYUV != true);
    INVALID_PARAM_IF(in_vi.Width < 16);
    INVALID_PARAM_IF(in_vi.Height < 16);
    INVALID_PARAM_IF(in_vi.Format.SSW < 0 || in_vi.Format.SSW > 4);
    INVALID_PARAM_IF(in_vi.Format.SSH < 0 || in_vi.Format.SSH > 4);
    INVALID_PARAM_IF(in_vi.Frames <= 0);
    INVALID_PARAM_IF(in_vi.Format.BitsPerSample < 8 || in_vi.Format.BitsPerSample > INTERNAL_BIT_DEPTH);
    INVALID_PARAM_IF(in_vi.Format.IsInteger != true);

    if (ep.output_depth < 0)
      ep.output_depth = in_vi.Format.BitsPerSample;
    if (ep.output_depth == 16)
        // set to appropriate precision mode
        ep.dither_algo = DA_16BIT_INTERLEAVED;

    int threshold_upper_limit = 64 * 8 - 1;
    int dither_upper_limit = 4096;

    #define CHECK_PARAM(value, lower_bound, upper_bound) \
    do { if ((int)value < (int)lower_bound || (int)value > (int)upper_bound) { snprintf(error_msg, sizeof(error_msg), "Invalid parameter %s, must be between %d and %d", #value, lower_bound, upper_bound); throw error_msg; } } while(0)

    CHECK_PARAM(ep.range, 0, 31);
    CHECK_PARAM(ep.Y, 0, threshold_upper_limit);
    CHECK_PARAM(ep.Cb, 0, threshold_upper_limit);
    CHECK_PARAM(ep.Cr, 0, threshold_upper_limit);
    CHECK_PARAM(ep.grainY, 0, dither_upper_limit);
    CHECK_PARAM(ep.grainC, 0, dither_upper_limit);
    CHECK_PARAM(ep.sample_mode, 1, 4);
    CHECK_PARAM(ep.dither_algo, DA_HIGH_NO_DITHERING, (DA_COUNT - 1) );
    CHECK_PARAM(ep.random_algo_ref, 0, (RANDOM_ALGORITHM_COUNT - 1) );
    CHECK_PARAM(ep.random_algo_grain, 0, (RANDOM_ALGORITHM_COUNT - 1) );
    

    // now the internal bit depth is 16, 
    // scale parameters to be consistent with 14bit range in previous versions
    ep.Y <<= 2;
    ep.Cb <<= 2;
    ep.Cr <<= 2;
    ep.grainY <<= 2;
    ep.grainC <<= 2;

    out_vi = in_vi;
    out_vi.Format.BitsPerSample = ep.output_depth;
    out_vi.Format.BytesPerSample = ep.output_depth == 8 ? 1 : 2;

    try
    {
        engine = new f3kdb_core_t(in_vi, ep, opt);
    } catch (std::bad_alloc&) {
        throw "Memory allocation failed";
    }
  }

  DSFrame GetFrame(int n, std::unordered_map<int, DSFrame> in_frames) override
  {
    auto src = in_frames[n];
    auto dst = src.Create(out_vi);
    std::for_each_n(std::execution::par_unseq, reinterpret_cast<char*>(0), in_vi.Format.Planes, [&](char&idx) {
      int p = static_cast<int>(reinterpret_cast<intptr_t>(&idx));
      auto src_stride = src.StrideBytes[p];
      auto src_ptr = src.SrcPointers[p];
      auto dst_stride = dst.StrideBytes[p];
      auto dst_ptr = dst.DstPointers[p];

      engine->process_plane(n, p, dst_ptr, dst_stride, src_ptr, src_stride);
    });

    return dst;
  }

  DSVideoInfo GetOutputVI()
  {
    return out_vi;
  }

  ~F3KDB() = default;
};
