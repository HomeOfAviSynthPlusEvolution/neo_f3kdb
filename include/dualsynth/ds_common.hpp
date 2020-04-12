/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - Common header+.
 *
 */
#pragma once

#include <avisynth.h>
#include <VapourSynth.h>
#include <cstring>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "ds_format.hpp"
#include "ds_videoinfo.hpp"
#include "ds_frame.hpp"

typedef void (*register_vsfilter_proc)(VSRegisterFunction, VSPlugin*);
typedef void (*register_avsfilter_proc)(IScriptEnvironment* env);
std::vector<register_vsfilter_proc> RegisterVSFilters();
std::vector<register_avsfilter_proc> RegisterAVSFilters();

enum ParamType
{
  Clip, Integer, Float, Boolean, String
};

struct Param
{
  const char* Name;
  const ParamType Type;
  const bool IsArray {false};
  bool AVSEnabled {true};
  bool VSEnabled {true};
  const bool IsOptional {true};
};

struct InDelegator
{
  virtual void Read(const char* name, int& output) = 0;
  virtual void Read(const char* name, int64_t& output) = 0;
  virtual void Read(const char* name, float& output) = 0;
  virtual void Read(const char* name, double& output) = 0;
  virtual void Read(const char* name, bool& output) = 0;
  virtual void Read(const char* name, std::string& output) = 0;
  virtual void Read(const char* name, std::vector<int>& output) = 0;
  virtual void Read(const char* name, std::vector<int64_t>& output) = 0;
  virtual void Read(const char* name, std::vector<float>& output) = 0;
  virtual void Read(const char* name, std::vector<double>& output) = 0;
  virtual void Read(const char* name, std::vector<bool>& output) = 0;
  virtual void Read(const char* name, void*& output) = 0;
  virtual void Free(void*& clip) = 0;
};

struct FetchFrameFunctor
{
  virtual DSFrame operator()(int n) = 0;
  virtual ~FetchFrameFunctor() {}
};
