/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - AviSynth+.
 *
 */

#pragma once

namespace Plugin {
  extern const char* Description;
}

namespace AVSInterface
{
  struct AVSInDelegator final : InDelegator {
    const AVSValue _args;
    std::unordered_map<std::string, int> _params_index_map;
    int NameToIndex(const char* name) {
      std::string name_string(name);
      if (_params_index_map.find(name_string) == _params_index_map.end())
        throw "Unknown parameter during NameToIndex";
      return _params_index_map[name_string];
    }
    void Read(const char* name, int& output) override {
      output = _args[NameToIndex(name)].AsInt(output);
    }
    void Read(const char* name, int64_t& output) override {
      output = _args[NameToIndex(name)].AsInt(static_cast<int>(output));
    }
    void Read(const char* name, float& output) override {
      output = static_cast<float>(_args[NameToIndex(name)].AsFloat(output));
    }
    void Read(const char* name, double& output) override {
      auto _default = output;
      output = _args[NameToIndex(name)].AsFloat(NAN);
      if (std::isnan(output))
        output = _default;
    }
    void Read(const char* name, bool& output) override {
      output = _args[NameToIndex(name)].AsBool(output);
    }
    void Read(const char* name, std::string& output) override {
      const char * result = _args[NameToIndex(name)].AsString(output.c_str());
      if (result)
        output = result;
    }
    void Read(const char* name, void*& output) override {
      PClip* clip = new PClip(_args[NameToIndex(name)].AsClip());
      output = (void *)(clip);
    }
    void Read(const char* name, std::vector<int>& output) override {
      auto arg = _args[NameToIndex(name)];
      if (!arg.IsArray())
        throw "Argument is not array";
      auto size = arg.ArraySize();
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(arg[i].AsInt());
    }
    void Read(const char* name, std::vector<int64_t>& output) override {
      auto arg = _args[NameToIndex(name)];
      if (!arg.IsArray())
        throw "Argument is not array";
      auto size = arg.ArraySize();
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(arg[i].AsInt());
    }
    void Read(const char* name, std::vector<float>& output) override {
      auto arg = _args[NameToIndex(name)];
      if (!arg.IsArray())
        throw "Argument is not array";
      auto size = arg.ArraySize();
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(static_cast<float>(arg[i].AsFloat()));
    }
    void Read(const char* name, std::vector<double>& output) override {
      auto arg = _args[NameToIndex(name)];
      if (!arg.IsArray())
        throw "Argument is not array";
      auto size = arg.ArraySize();
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(arg[i].AsFloat());
    }
    void Read(const char* name, std::vector<bool>& output) override {
      auto arg = _args[NameToIndex(name)];
      if (!arg.IsArray())
        throw "Argument is not array";
      auto size = arg.ArraySize();
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(arg[i].AsBool());
    }
    void Free(void*& clip) override {
      PClip* c = (PClip *)(clip);
      delete c;
      clip = nullptr;
    }

    AVSInDelegator(const AVSValue args, std::vector<Param> params) : _args(args)
    {
      int idx = 0;
      for (auto &&param : params)
      {
        if (!param.AVSEnabled) continue;
        _params_index_map[param.Name] = idx++;
      }
    }
  };

  struct AVSFetchFrameFunctor final : FetchFrameFunctor {
    PClip _clip;
    VideoInfo _vi;
    IScriptEnvironment* _env;
    std::mutex fetch_frame_mutex;
    AVSFetchFrameFunctor(PClip clip, VideoInfo vi, IScriptEnvironment * env)
      : _clip(clip), _vi(vi), _env(env) {}
    DSFrame operator()(int n) override {
      std::lock_guard<std::mutex> guard(fetch_frame_mutex);
      auto frame = _clip->GetFrame(n, _env);
      return DSFrame(frame, _vi, _env);
    }
    ~AVSFetchFrameFunctor() override {}
  };

  template<typename FilterType>
  struct AVSWrapper : IClip
  {
    AVSValue _args;
    IScriptEnvironment* _env;
    FilterType data;
    PClip clip;
    VideoInfo vi;
    AVSFetchFrameFunctor* functor {nullptr};
    
    AVSWrapper(AVSValue args, IScriptEnvironment* env)
      : _args(args), _env(env) {}
    
    void Initialize()
    {
      auto input_vi = DSVideoInfo();
      if (_args[0].IsClip()) {
        clip = _args[0].AsClip();
        input_vi = DSVideoInfo(clip->GetVideoInfo());
        functor = new AVSFetchFrameFunctor(clip, clip->GetVideoInfo(), _env);
      }
      auto argument = AVSInDelegator(_args, data.Params());
      data.Initialize(&argument, input_vi, functor);
    }

    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment * env) override {
      std::unordered_map<int, DSFrame> in_frames;
      if (functor) {
        std::vector<int> requests = data.RequestReferenceFrames(n);
        for (auto &&i : requests) {
          auto frame = clip->GetFrame(i, env);
          in_frames[i] = DSFrame(frame, vi, env);
        }
      }
      else
        in_frames[n] = DSFrame(env);
      
      return data.GetFrame(n, in_frames).ToAVSFrame();
    }

    const VideoInfo& __stdcall GetVideoInfo() override {
      auto output_vi = data.GetOutputVI();
      vi = output_vi.ToAVSVI();
      return vi;
    }

    void __stdcall GetAudio(void* buf, int64_t start, int64_t count, IScriptEnvironment* env) override { if (clip) clip->GetAudio(buf, start, count, env); }
    bool __stdcall GetParity(int n) override { return clip ? clip->GetParity(n) : false; }
    int __stdcall SetCacheHints(int cachehints, int frame_range) override { return data.SetCacheHints(cachehints, frame_range); }
    ~AVSWrapper() {
      delete functor;
    }
  };

  template<typename FilterType>
  AVSValue __cdecl Create(AVSValue args, void* user_data, IScriptEnvironment* env)
  {
    auto filter = new AVSWrapper<FilterType>(args, env);
    try {
      filter->Initialize();
    }
    catch (const char *err) {
      env->ThrowError("%s: %s", filter->data.AVSName(), err);
    }
    return filter;
  }

  template<typename FilterType>
  void RegisterFilter(IScriptEnvironment* env) {
    FilterType filter;
    env->AddFunction(filter.AVSName(), filter.AVSParams().c_str(), Create<FilterType>, nullptr);
  }
}

const AVS_Linkage *AVS_linkage = NULL;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, AVS_Linkage* linkage)
{
  AVS_linkage = linkage;
  auto filters = RegisterAVSFilters();
  for (auto &&RegisterFilter : filters) {
    RegisterFilter(env);
  }
  return Plugin::Description;
}
