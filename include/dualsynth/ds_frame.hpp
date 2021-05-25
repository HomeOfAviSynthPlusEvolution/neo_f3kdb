/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - DSFrame.
 *
 */

#pragma once

struct DSFrame
{
  int FrameWidth {0}, FrameHeight {0};

  const unsigned char ** SrcPointers {nullptr};
  int * StrideBytes {nullptr};
  unsigned char ** DstPointers {nullptr};
  DSFormat Format;

  // VapourSynth Interface
  const VSFrameRef* _vssrc {nullptr};
  VSFrameRef* _vsdst {nullptr};
  const VSCore* _vscore {nullptr};
  const VSAPI* _vsapi {nullptr};
  const VSFormat* _vsformat {nullptr};

  // AviSynth+ Interface
  PVideoFrame _avssrc;
  VideoInfo _vi;
  IScriptEnvironment * _env {nullptr};
  int planes_y[4] = { PLANAR_Y, PLANAR_U, PLANAR_V, PLANAR_A };
  int planes_r[4] = { PLANAR_R, PLANAR_G, PLANAR_B, PLANAR_A };
  int *planes {0};

  DSFrame() {}
  DSFrame(const VSCore* vscore, const VSAPI* vsapi)
    : _vscore(vscore), _vsapi(vsapi) {}
  DSFrame(const VSFrameRef* src, const VSCore* vscore, const VSAPI* vsapi)
    : _vssrc(src), _vscore(vscore), _vsapi(vsapi)
    , _vsformat(src ? _vsapi->getFrameFormat(src) : nullptr)
  {
    if (_vssrc) {
      Format = DSFormat(_vsformat);
      FrameWidth = _vsapi->getFrameWidth(src, 0);
      FrameHeight = _vsapi->getFrameHeight(src, 0);

      SrcPointers = new const unsigned char*[Format.Planes];
      StrideBytes = new int[Format.Planes];
      for (int i = 0; i < Format.Planes; i++) {
        SrcPointers[i] = _vsapi->getReadPtr(src, i);
        StrideBytes[i] = _vsapi->getStride(src, i);
      }
    }
  }

  DSFrame(IScriptEnvironment * env)
    : _env(env) {}
  DSFrame(PVideoFrame &src, VideoInfo vi, IScriptEnvironment * env)
    : _avssrc(src), _vi(vi), _env(env)
  {
    if (_avssrc) {
      Format = DSFormat(_vi.pixel_type);
      planes = Format.IsFamilyYUV ? planes_y : planes_r;
      FrameWidth = _vi.width;
      FrameHeight = _vi.height;

      SrcPointers = new const unsigned char*[Format.Planes];
      StrideBytes = new int[Format.Planes];
      for (int i = 0; i < Format.Planes; i++) {
        SrcPointers[i] = src->GetReadPtr(planes[i]);
        StrideBytes[i] = src->GetPitch(planes[i]);
      }
    }
  }

  DSFrame Create() { return Create(false, false); }
  DSFrame Create(bool copy) { return Create(copy, false); }
  DSFrame Create(bool copy, bool inplace)
  {
    if (_vssrc) {
      // Create a new VS frame
      const VSFrameRef* copy_frames[1] {ToVSFrame()};
      int copy_planes[4] = {0};
      auto vsframe = copy ?
        _vsapi->newVideoFrame2(_vsformat, FrameWidth, FrameHeight, copy_frames, copy_planes, ToVSFrame(), const_cast<VSCore*>(_vscore)) :
        _vsapi->newVideoFrame(_vsformat, FrameWidth, FrameHeight, ToVSFrame(), const_cast<VSCore*>(_vscore));
      _vsapi->freeFrame(copy_frames[0]);

      DSFrame new_frame(vsframe, _vscore, _vsapi);
      new_frame._vsdst = vsframe;
      new_frame.DstPointers = new unsigned char*[Format.Planes];
      for (int i = 0; i < Format.Planes; i++)
        new_frame.DstPointers[i] = _vsapi->getWritePtr(vsframe, i);
      return new_frame;
    }
    else if(_avssrc) {
      // Create a new AVS frame
      return Create(_vi);
    }
    throw "Unable to create from nothing.";
  }
  DSFrame Create(DSVideoInfo vi) {
    planes = vi.Format.IsFamilyYUV ? planes_y : planes_r;
    if (_vsapi) {
      auto vsframe = _vsapi->newVideoFrame(vi.Format.ToVSFormat(_vscore, _vsapi), vi.Width, vi.Height, ToVSFrame(), const_cast<VSCore*>(_vscore));
      DSFrame new_frame(vsframe, _vscore, _vsapi);
      new_frame._vsdst = vsframe;
      new_frame.DstPointers = new unsigned char*[Format.Planes];
      for (int i = 0; i < Format.Planes; i++)
        new_frame.DstPointers[i] = _vsapi->getWritePtr(vsframe, i);
      return new_frame;
    }
    else if (_env) {
      auto avsvi = vi.ToAVSVI();
      bool has_at_least_v8 = true;
      try { _env->CheckVersion(8); }
      catch (const AvisynthError&) { has_at_least_v8 = false; }
      auto new_avsframe = (has_at_least_v8) ? _env->NewVideoFrameP(avsvi, &_avssrc) : _env->NewVideoFrame(avsvi);
      auto dstp = new unsigned char*[Format.Planes];
      for (int i = 0; i < Format.Planes; i++)
        dstp[i] = new_avsframe->GetWritePtr(planes[i]);
      DSFrame new_frame(new_avsframe, avsvi, _env);
      new_frame.DstPointers = dstp;
      return new_frame;
    }
    throw "Unable to create from nothing.";
  }

  const VSFrameRef* ToVSFrame()
  {
    return _vsdst ? _vsapi->cloneFrameRef(_vsdst) :
           _vssrc ? _vsapi->cloneFrameRef(_vssrc) :
           nullptr;
  }
  PVideoFrame ToAVSFrame() {return _avssrc ? _avssrc : nullptr;}

  ~DSFrame()
  {
    if (SrcPointers)
      delete[] SrcPointers;
    if (DstPointers)
      delete[] DstPointers;
    if (StrideBytes)
      delete[] StrideBytes;
    if (_vsdst && _vsdst != _vssrc)
      _vsapi->freeFrame(_vsdst);
    if (_vssrc)
      _vsapi->freeFrame(_vssrc);
  }

  DSFrame(const DSFrame & old)
  {
    _avssrc = old._avssrc;
    std::memcpy(this, &old, sizeof(DSFrame));
    if (old.SrcPointers) {
      SrcPointers = new const unsigned char*[Format.Planes];
      std::copy_n(old.SrcPointers, Format.Planes, SrcPointers);
    }
    if (old.DstPointers) {
      DstPointers = new unsigned char*[Format.Planes];
      std::copy_n(old.DstPointers, Format.Planes, DstPointers);
    }
    if (old.StrideBytes) {
      StrideBytes = new int[Format.Planes];
      std::copy_n(old.StrideBytes, Format.Planes, StrideBytes);
    }
    if (_vsdst && _vsdst != _vssrc)
      _vsdst = const_cast<VSFrameRef*>(_vsapi->cloneFrameRef(old._vsdst));
    if (_vssrc)
      _vssrc = _vsapi->cloneFrameRef(old._vssrc);
  }
  DSFrame& operator =(const DSFrame & old)
  {
    if (&old == this)
      return *this;

    if (SrcPointers)
      delete[] SrcPointers;
    if (DstPointers)
      delete[] DstPointers;
    if (StrideBytes)
      delete[] StrideBytes;
    if (_vsdst && _vsdst != _vssrc)
      _vsapi->freeFrame(_vsdst);
    if (_vssrc)
      _vsapi->freeFrame(_vssrc);

    _avssrc = old._avssrc;
    std::memcpy(this, &old, sizeof(DSFrame));
    if (old.SrcPointers) {
      SrcPointers = new const unsigned char*[Format.Planes];
      std::copy_n(old.SrcPointers, Format.Planes, SrcPointers);
    }
    if (old.DstPointers) {
      DstPointers = new unsigned char*[Format.Planes];
      std::copy_n(old.DstPointers, Format.Planes, DstPointers);
    }
    if (old.StrideBytes) {
      StrideBytes = new int[Format.Planes];
      std::copy_n(old.StrideBytes, Format.Planes, StrideBytes);
    }
    if (_vsdst && _vsdst != _vssrc)
      _vsdst = const_cast<VSFrameRef*>(_vsapi->cloneFrameRef(old._vsdst));
    if (_vssrc)
      _vssrc = _vsapi->cloneFrameRef(old._vssrc);
    return *this;
  }
  DSFrame(DSFrame && old) noexcept
  {
    _avssrc = old._avssrc;
    std::memcpy(this, &old, sizeof(DSFrame));
    old.SrcPointers = nullptr;
    old.DstPointers = nullptr;
    old.StrideBytes = nullptr;
    old._vssrc = nullptr;
    old._vsdst = nullptr;
  }
  DSFrame& operator =(DSFrame && old) noexcept
  {
    if (&old == this)
      return *this;

    if (SrcPointers)
      delete[] SrcPointers;
    if (DstPointers)
      delete[] DstPointers;
    if (StrideBytes)
      delete[] StrideBytes;
    if (_vsdst && _vsdst != _vssrc)
      _vsapi->freeFrame(_vsdst);
    if (_vssrc)
      _vsapi->freeFrame(_vssrc);

    _avssrc = old._avssrc;
    std::memcpy(this, &old, sizeof(DSFrame));
    old.SrcPointers = nullptr;
    old.DstPointers = nullptr;
    old.StrideBytes = nullptr;
    old._vssrc = nullptr;
    old._vsdst = nullptr;
    return *this;
  }
};
