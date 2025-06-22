/*
 * Copyright 2020 Xinyue Lu
 *
 * Temporal Median - filter.
 *
 */

#pragma once

#include <memory>

#ifdef HAS_EXECUTION
  #include <execution>
#endif

#ifndef __cpp_lib_execution
  #undef ENABLE_PAR
#endif

#ifdef ENABLE_PAR
  #define PAR_POLICY std::execution::par
#else
  #define PAR_POLICY nullptr
#endif

#include "compiler_compat.h"
#include "core.h"
#include "constants.h"
#include "impl_dispatch.h"

int GetCPUFlags();

struct F3KDB final : Filter {
  f3kdb_params_t ep;
  std::unique_ptr<f3kdb_core_t> engine;
  char error_msg[1024];
  DSVideoInfo out_vi;
  bool mt {true};

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
      Param {"opt", Integer},
      Param {"mt", Boolean},
      Param {"dither_algo", Integer},
      Param {"keep_tv_range", Boolean},
      Param {"output_depth", Integer},
      Param {"random_algo_ref", Integer},
      Param {"random_algo_grain", Integer},
      Param {"random_param_ref", Float},
      Param {"random_param_grain", Float},
      Param {"preset", String},
      Param{ "y_1", Integer},
      Param{ "cb_1", Integer},
      Param{ "cr_1", Integer},
      Param{ "y_2", Integer },
      Param{ "cb_2", Integer },
      Param{ "cr_2", Integer },
      Param{ "scale", Boolean },
      Param{ "angle_boost", Float },
      Param{ "max_angle", Float },
    };
  }
  void Initialize(InDelegator* in, DSVideoInfo in_vi, FetchFrameFunctor* fetch_frame) override
  {
    Filter::Initialize(in, in_vi, fetch_frame);
    std::string preset;
    in->Read("preset", preset);
    std::istringstream piss(preset);

    bool scale = false;
    in->Read("scale", scale);

    while(!piss.eof()) {
      std::string piss1;
      std::getline(piss, piss1, '/');
      if (piss1 == "depth")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = 0;
      else if (piss1 == "low")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = (scale) ? 128 : 32;
      else if (piss1 == "medium")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = (scale) ? 192 : 48;
      else if (piss1 == "high")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = (scale) ? 256 : 64;
      else if (piss1 == "veryhigh")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = (scale) ? 320 : 80;
      else if (piss1 == "nograin")
        ep.grainY = ep.grainC = 0;
      else if (piss1 == "luma")
        ep.Cb = ep.Cr = ep.grainC = 0;
      else if (piss1 == "chroma")
        ep.Y = ep.grainY = 0;
    }
    int tmp;
    in->Read("range", ep.range);
    in->Read("y", ep.Y);
    in->Read("cb", ep.Cb);
    in->Read("cr", ep.Cr);
    in->Read("grainy", ep.grainY);
    in->Read("grainc", ep.grainC);
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
    in->Read("y_1", ep.Y_1);
    in->Read("cb_1", ep.Cb_1);
    in->Read("cr_1", ep.Cr_1);
    in->Read("y_2", ep.Y_2);
    in->Read("cb_2", ep.Cb_2);
    in->Read("cr_2", ep.Cr_2);
    in->Read("angle_boost", ep.angle_boost);
    in->Read("max_angle", ep.max_angle);

    ep.Y_1 = ep.Y_1 == -1 ? ep.Y : ep.Y_1;
    ep.Cb_1 = ep.Cb_1 == -1 ? ep.Cb : ep.Cb_1;
    ep.Cr_1 = ep.Cr_1 == -1 ? ep.Cr : ep.Cr_1;
    ep.Y_2 = ep.Y_2 == -1 ? ep.Y : ep.Y_2;
    ep.Cb_2 = ep.Cb_2 == -1 ? ep.Cb : ep.Cb_2;
    ep.Cr_2 = ep.Cr_2 == -1 ? ep.Cr : ep.Cr_2;

    int opt_in = -1;
    in->Read("opt", opt_in);
    in->Read("mt", mt);

    OPTIMIZATION_MODE opt = [&]() {
        const int CPUFlags = GetCPUFlags();

        if (ep.sample_mode >= 5 && ep.sample_mode <= 7) {
            const int AVX512_REQUIRED_FLAGS = CPUF_AVX512F | CPUF_AVX512BW | CPUF_AVX512DQ | CPUF_AVX512VL | CPUF_AVX512CD;

            if (((CPUFlags & AVX512_REQUIRED_FLAGS) == AVX512_REQUIRED_FLAGS) && (opt_in == 3 || opt_in < 0))
                return IMPL_AVX512;

            if ((CPUFlags & CPUF_AVX2) && (opt_in == 2 || opt_in < 0))
                return IMPL_AVX2;
        }

        if ((CPUFlags & CPUF_SSE4_1) && (opt_in > 0 || opt_in < 0))
            return IMPL_SSE4;

        return IMPL_C;
        }();

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

    const int y_threshold_upper_limit = scale ? 65535 : 511;
    const int cb_threshold_upper_limit = scale ? 65535 : 511;
    const int cr_threshold_upper_limit = scale ? 65535 : 511;
    constexpr int dither_upper_limit = 4096;

    #define CHECK_PARAM(value, lower_bound, upper_bound) \
    do { if ((int)value < (int)lower_bound || (int)value > (int)upper_bound) { snprintf(error_msg, sizeof(error_msg), "Invalid parameter %s, must be between %d and %d", #value, lower_bound, upper_bound); throw error_msg; } } while(0)

    CHECK_PARAM(ep.range, 0, 255);
    CHECK_PARAM(ep.Y, 0, y_threshold_upper_limit);
    CHECK_PARAM(ep.Cb, 0, cb_threshold_upper_limit);
    CHECK_PARAM(ep.Cr, 0, cr_threshold_upper_limit);
    CHECK_PARAM(ep.grainY, 0, dither_upper_limit);
    CHECK_PARAM(ep.grainC, 0, dither_upper_limit);
    CHECK_PARAM(ep.sample_mode, 1, 7);
    CHECK_PARAM(ep.dither_algo, DA_HIGH_NO_DITHERING, (DA_COUNT - 1) );
    CHECK_PARAM(ep.random_algo_ref, 0, (RANDOM_ALGORITHM_COUNT - 1) );
    CHECK_PARAM(ep.random_algo_grain, 0, (RANDOM_ALGORITHM_COUNT - 1) );
    CHECK_PARAM(ep.Y_1, 0, y_threshold_upper_limit);
    CHECK_PARAM(ep.Cb_1, 0, cb_threshold_upper_limit);
    CHECK_PARAM(ep.Cr_1, 0, cr_threshold_upper_limit);
    CHECK_PARAM(ep.Y_2, 0, y_threshold_upper_limit);
    CHECK_PARAM(ep.Cb_2, 0, cb_threshold_upper_limit);
    CHECK_PARAM(ep.Cr_2, 0, cr_threshold_upper_limit);

    if (ep.angle_boost < 0.0f)
        throw "invalid parameter angle_boost, must be positive value";

    if (ep.max_angle < 0.0f || ep.max_angle > 1.0f)
        throw "invalid parameter max_angle, must be between 0.0 and 1.0";

    // now the internal bit depth is 16,
    // scale parameters to be consistent with 14bit range in previous versions
    ep.Y = scale ? ep.Y : ep.Y << 2;
    ep.Cb = scale ? ep.Cb : ep.Cb << 2;
    ep.Cr = scale ? ep.Cr : ep.Cr << 2;
    ep.Y_1 = scale ? ep.Y_1 : ep.Y_1 << 2;
    ep.Cb_1 = scale ? ep.Cb_1 : ep.Cb_1 << 2;
    ep.Cr_1 = scale ? ep.Cr_1 : ep.Cr_1 << 2;
    ep.Y_2 = scale ? ep.Y_2 : ep.Y_2 << 2;
    ep.Cb_2 = scale ? ep.Cb_2 : ep.Cb_2 << 2;
    ep.Cr_2 = scale ? ep.Cr_2 : ep.Cr_2 << 2;
    ep.grainY <<= 2;
    ep.grainC <<= 2;

    out_vi = in_vi;
    out_vi.Format.BitsPerSample = ep.output_depth;
    out_vi.Format.BytesPerSample = ep.output_depth == 8 ? 1 : 2;

    try
    {
        engine = std::make_unique<f3kdb_core_t>(in_vi, ep, opt);
    } catch (std::bad_alloc&) {
        throw "Memory allocation failed";
    }
  }

  DSFrame GetFrame(int n, std::unordered_map<int, DSFrame> in_frames) override
  {
    auto src = in_frames[n];
    auto dst = src.Create(out_vi);
    auto core = [&](char&idx) {
      int p = static_cast<int>(reinterpret_cast<intptr_t>(&idx));
      auto src_stride = src.StrideBytes[p];
      auto src_ptr = src.SrcPointers[p];
      auto dst_stride = dst.StrideBytes[p];
      auto dst_ptr = dst.DstPointers[p];

      engine->process_plane(n, p, dst_ptr, dst_stride, src_ptr, src_stride);
    };

#ifdef ENABLE_PAR
    if(mt)
      std::for_each_n(PAR_POLICY, reinterpret_cast<char*>(0), in_vi.Format.Planes, core);
    else
#endif
    for (intptr_t i = 0; i < in_vi.Format.Planes; i++)
      core(*reinterpret_cast<char*>(i));

    return dst;
  }

  DSVideoInfo GetOutputVI() override
  {
    return out_vi;
  }

  ~F3KDB() = default;
};
