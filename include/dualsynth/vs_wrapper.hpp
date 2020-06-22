/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - VapourSynth.
 *
 */

#pragma once

namespace Plugin {
  extern const char* Identifier;
  extern const char* Namespace;
  extern const char* Description;
}

namespace VSInterface {
  const VSAPI * API;

  struct VSInDelegator final : InDelegator {
    const VSMap *_in;
    const VSAPI *_vsapi;
    int _err;
    void Read(const char* name, int& output) override {
      auto _default = output;
      output = static_cast<int>(_vsapi->propGetInt(_in, name, 0, &_err));
      if (_err) output = _default;
    }
    void Read(const char* name, int64_t& output) override {
      auto _default = output;
      output = _vsapi->propGetInt(_in, name, 0, &_err);
      if (_err) output = _default;
    }
    void Read(const char* name, float& output) override {
      auto _default = output;
      output = static_cast<float>(_vsapi->propGetFloat(_in, name, 0, &_err));
      if (_err) output = _default;
    }
    void Read(const char* name, double& output) override {
      auto _default = output;
      output = _vsapi->propGetFloat(_in, name, 0, &_err);
      if (_err) output = _default;
    }
    void Read(const char* name, bool& output) override {
      auto output_int = _vsapi->propGetInt(_in, name, 0, &_err);
      if (!_err) output = output_int != 0;
    }
    void Read(const char* name, std::string& output) override {
      auto output_str = _vsapi->propGetData(_in, name, 0, &_err);
      if (!_err) output = output_str;
    }
    void Read(const char* name, std::vector<int>& output) override {
      auto size = _vsapi->propNumElements(_in, name);
      if (size < 0) return;
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(static_cast<int>(_vsapi->propGetInt(_in, name, i, &_err)));
    }
    void Read(const char* name, std::vector<int64_t>& output) override {
      auto size = _vsapi->propNumElements(_in, name);
      if (size < 0) return;
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(_vsapi->propGetInt(_in, name, i, &_err));
    }
    void Read(const char* name, std::vector<float>& output) override {
      auto size = _vsapi->propNumElements(_in, name);
      if (size < 0) return;
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(static_cast<float>(_vsapi->propGetFloat(_in, name, i, &_err)));
    }
    void Read(const char* name, std::vector<double>& output) override {
      auto size = _vsapi->propNumElements(_in, name);
      if (size < 0) return;
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(_vsapi->propGetFloat(_in, name, i, &_err));
    }
    void Read(const char* name, std::vector<bool>& output) override {
      auto size = _vsapi->propNumElements(_in, name);
      if (size < 0) return;
      output.clear();
      for (int i = 0; i < size; i++)
        output.push_back(_vsapi->propGetInt(_in, name, i, &_err));
    }
    void Read(const char* name, void*& output) override {
      output = reinterpret_cast<void *>(_vsapi->propGetNode(_in, name, 0, &_err));
    }
    void Free(void*& clip) override {
      _vsapi->freeNode(reinterpret_cast<VSNodeRef *>(clip));
      clip = nullptr;
    }
    VSInDelegator(const VSMap *in, const VSAPI *vsapi) : _in(in), _vsapi(vsapi) {}
  };

  struct VSFetchFrameFunctor final : FetchFrameFunctor {
    VSNodeRef *_vs_clip;
    VSCore *_core;
    const VSAPI *_vsapi;
    VSFrameContext *_frameCtx;
    VSFetchFrameFunctor(VSNodeRef *clip, VSCore *core, const VSAPI *vsapi)
      : _vs_clip(clip), _core(core), _vsapi(vsapi) {}
    DSFrame operator()(int n) override {
      return DSFrame(_vsapi->getFrameFilter(n, _vs_clip, _frameCtx), _core, _vsapi);
    }
    ~VSFetchFrameFunctor() override {
      _vsapi->freeNode(_vs_clip);
    }
  };

  template<typename FilterType>
  void VS_CC Initialize(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) {
    auto Data = reinterpret_cast<FilterType*>(*instanceData);
    auto output_vi = Data->GetOutputVI();
    vsapi->setVideoInfo(output_vi.ToVSVI(core, vsapi), 1, node);
  }

  template<typename FilterType>
  void VS_CC Delete(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    auto filter = reinterpret_cast<FilterType*>(instanceData);
    auto functor = reinterpret_cast<VSFetchFrameFunctor*>(filter->fetch_frame);
    delete functor;
    delete filter;
  }

  template<typename FilterType>
  const VSFrameRef* VS_CC GetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    auto filter = reinterpret_cast<FilterType*>(*instanceData);
    auto functor = reinterpret_cast<VSFetchFrameFunctor*>(filter->fetch_frame);
    if (functor)
      functor->_frameCtx = frameCtx;

    std::vector<int> ref_frames;
    if (activationReason == VSActivationReason::arInitial) {
      if (functor) {
        ref_frames = filter->RequestReferenceFrames(n);
        for (auto &&i : ref_frames)
          vsapi->requestFrameFilter(i, functor->_vs_clip, frameCtx);
      }
      else {
        std::unordered_map<int, DSFrame> in_frames;
        in_frames[n] = DSFrame(core, vsapi);
        auto vs_frame = (filter->GetFrame(n, in_frames).ToVSFrame());
        return vs_frame;
      }
    }
    else if (activationReason == VSActivationReason::arAllFramesReady) {
      std::unordered_map<int, DSFrame> in_frames;
      if (functor) {
        ref_frames = filter->RequestReferenceFrames(n);
        for (auto &&i : ref_frames)
          in_frames[i] = DSFrame(vsapi->getFrameFilter(i, functor->_vs_clip, frameCtx), core, vsapi);
      }
      else
        in_frames[n] = DSFrame(core, vsapi);

      auto vs_frame = (filter->GetFrame(n, in_frames).ToVSFrame());
      return vs_frame;
    }
    return nullptr;
  }

  template<typename FilterType>
  void VS_CC Create(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    auto filter = new FilterType{};
    auto argument = VSInDelegator(in, vsapi);
    try {
      void* clip = nullptr;
      VSFetchFrameFunctor* functor = nullptr;
      DSVideoInfo input_vi;
      try {
        argument.Read("clip", clip);
        if (clip) {
          auto vs_clip = reinterpret_cast<VSNodeRef*>(clip);
          functor = new VSFetchFrameFunctor(vs_clip, core, vsapi);
          input_vi = DSVideoInfo(vsapi->getVideoInfo(vs_clip));
        }
      }
      catch(const char *) { /* No clip, source filter */ }
      filter->Initialize(&argument, input_vi, functor);
      vsapi->createFilter(in, out, filter->VSName(), Initialize<FilterType>, GetFrame<FilterType>, Delete<FilterType>, filter->VSMode(), 0, filter, core);
    }
    catch(const char *err){
      char msg_buff[256];
      snprintf(msg_buff, 256, "%s: %s", filter->VSName(), err);
      vsapi->setError(out, msg_buff);
      delete filter;
    }
  }

  template<typename FilterType>
  void RegisterFilter(VSRegisterFunction registerFunc, VSPlugin* vsplugin) {
    FilterType filter;
    registerFunc(filter.VSName(), filter.VSParams().c_str(), Create<FilterType>, nullptr, vsplugin);
  }

  void RegisterPlugin(VSConfigPlugin configFunc, VSPlugin* vsplugin) {
    configFunc(Plugin::Identifier, Plugin::Namespace, Plugin::Description, VAPOURSYNTH_API_VERSION, 1, vsplugin);
  }
}

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin* vsplugin) {
  VSInterface::RegisterPlugin(configFunc, vsplugin);
  auto filters = RegisterVSFilters();
  for (auto &&RegisterFilter : filters) {
    RegisterFilter(registerFunc, vsplugin);
  }
}
