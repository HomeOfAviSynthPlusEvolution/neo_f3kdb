/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - DSFrame.
 *
 */

#pragma once

#include "ds_format.hpp"
#include "ds_videoinfo.hpp"
#include <avisynth.h>
#include <VapourSynth4.h>
#include <cstring>
#include <algorithm>

struct DSFrame
{
	int FrameWidth{ 0 };
	int FrameHeight{ 0 };

	const uint8_t** SrcPointers{ nullptr };
	ptrdiff_t* StrideBytes{ nullptr };
	uint8_t** DstPointers{ nullptr };
	DSFormat Format;

	const VSFrame* _vssrc{ nullptr };
	VSFrame* _vsdst{ nullptr };
	VSCore* _vscore{ nullptr };
	const VSAPI* _vsapi{ nullptr };

	PVideoFrame _avssrc;
	VideoInfo _vi;
	IScriptEnvironment* _env{ nullptr };
	int planes_y[4] = { PLANAR_Y, PLANAR_U, PLANAR_V, PLANAR_A };
	int planes_r[4] = { PLANAR_R, PLANAR_G, PLANAR_B, PLANAR_A };
	int* planes{ planes_y };

	DSFrame() = default;

	DSFrame(VSCore* core, const VSAPI* vsapi)
		: _vscore(core), _vsapi(vsapi)
	{
	}

	DSFrame(const VSFrame* src, VSCore* core, const VSAPI* vsapi)
		: _vssrc(src), _vscore(core), _vsapi(vsapi)
	{
		if (_vssrc && _vsapi) {
			const VSVideoFormat* vsFormat = _vsapi->getVideoFrameFormat(_vssrc);
			if (!vsFormat)
				throw "Invalid frame format in DSFrame constructor (not a video frame)";
			Format = DSFormat(*vsFormat);
			FrameWidth = _vsapi->getFrameWidth(_vssrc, 0);
			FrameHeight = _vsapi->getFrameHeight(_vssrc, 0);

			SrcPointers = new const uint8_t * [Format.Planes];
			StrideBytes = new ptrdiff_t[Format.Planes];
			for (int i = 0; i < Format.Planes; i++) {
				SrcPointers[i] = _vsapi->getReadPtr(_vssrc, i);
				StrideBytes[i] = _vsapi->getStride(_vssrc, i);
			}
			planes = (Format.IsFamilyYUV || Format.IsFamilyGray) ? planes_y : planes_r;
		}
	}

	DSFrame(IScriptEnvironment* env)
		: _env(env)
	{
	}

	DSFrame(const PVideoFrame& src, VideoInfo vi, IScriptEnvironment* env)
		: _avssrc(src), _vi(vi), _env(env)
	{
		if (_avssrc) {
			Format = DSFormat(_vi);
			planes = (Format.IsFamilyYUV || Format.IsFamilyGray) ? planes_y : planes_r;
			FrameWidth = _vi.width;
			FrameHeight = _vi.height;

			SrcPointers = new const uint8_t * [Format.Planes];
			StrideBytes = new ptrdiff_t[Format.Planes];
			for (int i = 0; i < Format.Planes; i++) {
				SrcPointers[i] = src->GetReadPtr(planes[i]);
				StrideBytes[i] = src->GetPitch(planes[i]);
			}
		}
	}

	DSFrame Create() { return Create(false, false); }
	DSFrame Create(bool copy) { return Create(copy, false); }

	DSFrame Create(bool copy, bool /* inplace */)
	{
		if (_vssrc && _vsapi) {
			const VSVideoFormat* srcFormat = _vsapi->getVideoFrameFormat(_vssrc);
			if (!srcFormat)
				throw "Invalid source video format in DSFrame::Create";

			VSFrame* vsframe = nullptr;

			if (copy) {
				const VSFrame* planeSrc[4];
				int planeIndices[4];
				int numPlanes = (srcFormat->numPlanes > 4) ? 4 : srcFormat->numPlanes;

				for (int i = 0; i < numPlanes; i++) {
					planeSrc[i] = _vssrc;
					planeIndices[i] = i;
				}

				vsframe = _vsapi->newVideoFrame2(srcFormat, FrameWidth, FrameHeight,
					planeSrc, planeIndices, _vssrc, _vscore);
			}
			else {
				vsframe = _vsapi->newVideoFrame(srcFormat, FrameWidth, FrameHeight, _vssrc, _vscore);
			}

			if (!vsframe)
				throw "Failed to create VS frame";

			DSFrame new_frame;
			new_frame._vsdst = vsframe;
			new_frame._vscore = _vscore;
			new_frame._vsapi = _vsapi;
			new_frame.Format = DSFormat(*srcFormat);
			new_frame.FrameWidth = FrameWidth;
			new_frame.FrameHeight = FrameHeight;
			new_frame.SrcPointers = new const uint8_t * [new_frame.Format.Planes];
			new_frame.DstPointers = new uint8_t * [new_frame.Format.Planes];
			new_frame.StrideBytes = new ptrdiff_t[new_frame.Format.Planes];
			for (int i = 0; i < new_frame.Format.Planes; i++) {
				new_frame.DstPointers[i] = _vsapi->getWritePtr(vsframe, i);
				new_frame.SrcPointers[i] = new_frame.DstPointers[i];
				new_frame.StrideBytes[i] = _vsapi->getStride(vsframe, i);
			}
			new_frame.planes = (new_frame.Format.IsFamilyYUV || new_frame.Format.IsFamilyGray) ? new_frame.planes_y : new_frame.planes_r;

			return new_frame;
		}
		else if (_avssrc) {
			DSFrame new_frame = Create(_vi);

			if (copy) {
				for (int p = 0; p < Format.Planes; ++p) {
					int row_size = _avssrc->GetRowSize(planes[p]);
					int height = _avssrc->GetHeight(planes[p]);
					ptrdiff_t src_pitch = StrideBytes[p];
					ptrdiff_t dst_pitch = new_frame.StrideBytes[p];
					const uint8_t* src_ptr = SrcPointers[p];
					uint8_t* dst_ptr = new_frame.DstPointers[p];

					if (src_pitch == dst_pitch && row_size == src_pitch)
						std::memcpy(dst_ptr, src_ptr, src_pitch * height);
					else {
						for (int y = 0; y < height; ++y) {
							std::memcpy(dst_ptr, src_ptr, row_size);
							src_ptr += src_pitch;
							dst_ptr += dst_pitch;
						}
					}
				}
			}
			return new_frame;
		}
		throw "Unable to create from nothing.";
	}

	DSFrame Create(DSVideoInfo vi)
	{
		planes = (vi.Format.IsFamilyYUV || vi.Format.IsFamilyGray) ? planes_y : planes_r;

		if (_vsapi) {
			VSVideoFormat vsFormat{};
			if (!vi.Format.ToVSFormat(&vsFormat, _vscore, _vsapi))
				throw "Unable to convert DSFormat to VSVideoFormat";

			VSFrame* vsframe = _vsapi->newVideoFrame(&vsFormat, vi.Width, vi.Height, _vssrc, _vscore);
			if (!vsframe)
				throw "Failed to create VS frame from DSVideoInfo";

			DSFrame new_frame;
			new_frame._vsdst = vsframe;
			new_frame._vscore = _vscore;
			new_frame._vsapi = _vsapi;
			new_frame.Format = vi.Format;
			new_frame.FrameWidth = vi.Width;
			new_frame.FrameHeight = vi.Height;
			new_frame.DstPointers = new uint8_t * [new_frame.Format.Planes];
			new_frame.StrideBytes = new ptrdiff_t[new_frame.Format.Planes];
			for (int i = 0; i < new_frame.Format.Planes; i++) {
				new_frame.DstPointers[i] = _vsapi->getWritePtr(vsframe, i);
				new_frame.StrideBytes[i] = _vsapi->getStride(vsframe, i);
			}
			new_frame.planes = (new_frame.Format.IsFamilyYUV || new_frame.Format.IsFamilyGray) ? new_frame.planes_y : new_frame.planes_r;

			return new_frame;
		}
		else if (_env) {
			auto avsvi = vi.ToAVSVI();
			bool has_at_least_v8 = true;
			try { _env->CheckVersion(8); }
			catch (const AvisynthError&) { has_at_least_v8 = false; }

			auto new_avsframe = (has_at_least_v8)
				? _env->NewVideoFrameP(avsvi, _avssrc ? &_avssrc : nullptr)
				: _env->NewVideoFrame(avsvi);

			auto dstp = new uint8_t * [vi.Format.Planes];
			for (int i = 0; i < vi.Format.Planes; i++)
				dstp[i] = new_avsframe->GetWritePtr(planes[i]);

			DSFrame new_frame(new_avsframe, avsvi, _env);
			new_frame.DstPointers = dstp;
			return new_frame;
		}
		throw "Unable to create from nothing.";
	}

	const VSFrame* ToVSFrame() const
	{
		if (_vsapi) {
			if (_vsdst)
				return _vsapi->addFrameRef(_vsdst);
			if (_vssrc)
				return _vsapi->addFrameRef(_vssrc);
		}
		return nullptr;
	}

	PVideoFrame ToAVSFrame() const
	{
		return _avssrc;
	}

	~DSFrame()
	{
		if (SrcPointers) {
			delete[] SrcPointers;
			SrcPointers = nullptr;
		}
		if (DstPointers) {
			delete[] DstPointers;
			DstPointers = nullptr;
		}
		if (StrideBytes) {
			delete[] StrideBytes;
			StrideBytes = nullptr;
		}

		if (_vsapi) {
			if (_vsdst && _vsdst != _vssrc) {
				_vsapi->freeFrame(_vsdst);
				_vsdst = nullptr;
			}
			if (_vssrc) {
				_vsapi->freeFrame(_vssrc);
				_vssrc = nullptr;
			}
		}
	}

	DSFrame(const DSFrame& old)
		: FrameWidth(old.FrameWidth),
		FrameHeight(old.FrameHeight),
		Format(old.Format),
		_vscore(old._vscore),
		_vsapi(old._vsapi),
		_avssrc(old._avssrc),
		_vi(old._vi),
		_env(old._env)
	{
		if (_vsapi) {
			_vssrc = old._vssrc ? _vsapi->addFrameRef(old._vssrc) : nullptr;
			_vsdst = old._vsdst ? const_cast<VSFrame*>(_vsapi->addFrameRef(old._vsdst)) : nullptr;
		}
		else {
			_vssrc = nullptr;
			_vsdst = nullptr;
		}

		if (old.SrcPointers) {
			SrcPointers = new const uint8_t * [Format.Planes];
			std::copy_n(old.SrcPointers, Format.Planes, SrcPointers);
		}
		if (old.DstPointers) {
			DstPointers = new uint8_t * [Format.Planes];
			std::copy_n(old.DstPointers, Format.Planes, DstPointers);
		}
		if (old.StrideBytes) {
			StrideBytes = new ptrdiff_t[Format.Planes];
			std::copy_n(old.StrideBytes, Format.Planes, StrideBytes);
		}

		planes = (Format.IsFamilyYUV || Format.IsFamilyGray) ? planes_y : planes_r;
	}

	DSFrame& operator=(const DSFrame& old)
	{
		if (&old == this) return *this;

		delete[] SrcPointers; SrcPointers = nullptr;
		delete[] DstPointers; DstPointers = nullptr;
		delete[] StrideBytes; StrideBytes = nullptr;

		if (_vsapi) {
			if (_vsdst && _vsdst != _vssrc) _vsapi->freeFrame(_vsdst);
			if (_vssrc) _vsapi->freeFrame(_vssrc);
		}

		_avssrc = old._avssrc;
		_vi = old._vi;
		_env = old._env;
		Format = old.Format;
		FrameWidth = old.FrameWidth;
		FrameHeight = old.FrameHeight;
		_vscore = old._vscore;
		_vsapi = old._vsapi;

		if (_vsapi) {
			_vssrc = old._vssrc ? _vsapi->addFrameRef(old._vssrc) : nullptr;
			_vsdst = old._vsdst ? const_cast<VSFrame*>(_vsapi->addFrameRef(old._vsdst)) : nullptr;
		}
		else {
			_vssrc = nullptr;
			_vsdst = nullptr;
		}

		if (old.SrcPointers) {
			SrcPointers = new const uint8_t * [Format.Planes];
			std::copy_n(old.SrcPointers, Format.Planes, SrcPointers);
		}
		if (old.DstPointers) {
			DstPointers = new uint8_t * [Format.Planes];
			std::copy_n(old.DstPointers, Format.Planes, DstPointers);
		}
		if (old.StrideBytes) {
			StrideBytes = new ptrdiff_t[Format.Planes];
			std::copy_n(old.StrideBytes, Format.Planes, StrideBytes);
		}

		planes = (Format.IsFamilyYUV || Format.IsFamilyGray) ? planes_y : planes_r;
		return *this;
	}

	DSFrame(DSFrame&& old) noexcept
		: FrameWidth(old.FrameWidth),
		FrameHeight(old.FrameHeight),
		SrcPointers(old.SrcPointers),
		StrideBytes(old.StrideBytes),
		DstPointers(old.DstPointers),
		Format(old.Format),
		_vssrc(old._vssrc),
		_vsdst(old._vsdst),
		_vscore(old._vscore),
		_vsapi(old._vsapi),
		_avssrc(std::move(old._avssrc)),
		_vi(old._vi),
		_env(old._env)
	{
		old.SrcPointers = nullptr;
		old.DstPointers = nullptr;
		old.StrideBytes = nullptr;
		old._vssrc = nullptr;
		old._vsdst = nullptr;

		planes = (Format.IsFamilyYUV || Format.IsFamilyGray) ? planes_y : planes_r;
	}

	DSFrame& operator=(DSFrame&& old) noexcept
	{
		if (&old == this) return *this;

		delete[] SrcPointers;
		delete[] DstPointers;
		delete[] StrideBytes;

		if (_vsapi) {
			if (_vsdst && _vsdst != _vssrc) _vsapi->freeFrame(_vsdst);
			if (_vssrc) _vsapi->freeFrame(_vssrc);
		}

		_avssrc = std::move(old._avssrc);
		_vi = old._vi;
		_env = old._env;
		Format = old.Format;
		FrameWidth = old.FrameWidth;
		FrameHeight = old.FrameHeight;
		_vscore = old._vscore;
		_vsapi = old._vsapi;
		_vssrc = old._vssrc;
		_vsdst = old._vsdst;
		SrcPointers = old.SrcPointers;
		DstPointers = old.DstPointers;
		StrideBytes = old.StrideBytes;

		old.SrcPointers = nullptr;
		old.DstPointers = nullptr;
		old.StrideBytes = nullptr;
		old._vssrc = nullptr;
		old._vsdst = nullptr;
		old._avssrc = PVideoFrame();

		planes = (Format.IsFamilyYUV || Format.IsFamilyGray) ? planes_y : planes_r;
		return *this;
	}
};
