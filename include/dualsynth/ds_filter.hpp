/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - Filter parent class.
 *
 */

#pragma once

struct Filter
{
  DSVideoInfo in_vi;
  FetchFrameFunctor* fetch_frame;
  virtual const char* VSName() const { return "FilterFoo"; }
  virtual const char* AVSName() const { return "FilterFoo"; }
  virtual const MtMode AVSMode() const { return MT_SERIALIZED; }
  virtual const VSFilterMode VSMode() const { return fmSerial; }
  virtual const std::vector<Param> Params() const = 0;
  virtual const std::string VSParams() const
  {
    std::stringstream ss;
    auto params = this->Params();
    for (auto &&p : params)
    {
      if (!p.VSEnabled) continue;
      std::string type_name;
      switch(p.Type) {
        case Clip: type_name = "clip"; break;
        case Integer: type_name = "int"; break;
        case Float: type_name = "float"; break;
        case Boolean: type_name = "int"; break;
        case String: type_name = "data"; break;
      }
      ss << p.Name << ':' << type_name;
      if (p.IsArray)
        ss << "[]";
      if (p.IsOptional)
        ss << ":opt";
      ss << ';';
    }
    return ss.str();
  };
  virtual const std::string AVSParams() const
  {
    std::stringstream ss;
    auto params = this->Params();
    for (auto &&p : params)
    {
      if (!p.AVSEnabled) continue;
      char type_name;
      switch(p.Type) {
        case Clip: type_name = 'c'; break;
        case Integer: type_name = 'i'; break;
        case Float: type_name = 'f'; break;
        case Boolean: type_name = 'b'; break;
        case String: type_name = 's'; break;
      }
      if (p.IsOptional)
        ss << '[' << p.Name << ']';
      ss << type_name;
    }
    return ss.str();
  };
  virtual void Initialize(InDelegator* in, DSVideoInfo in_vi, FetchFrameFunctor* fetch_frame)
  {
    this->in_vi = in_vi;
    this->fetch_frame = fetch_frame;
  };
  virtual std::vector<int> RequestReferenceFrames(int n) const
  {
    return std::vector<int>{n};
  }
  virtual DSFrame GetFrame(int n, std::unordered_map<int, DSFrame> in_frames)
  {
    return in_frames.size() > 0 ? in_frames.begin()->second : DSFrame();
  }
  virtual DSVideoInfo GetOutputVI()
  {
    return in_vi;
  }
  virtual int SetCacheHints(int cachehints, int frame_range)
  {
    return cachehints == CACHE_GET_MTMODE ? AVSMode() : 0;
  }
};
