/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - DSVideoInfo.
 *
 */

#pragma once

#include "ds_format.hpp"
#include <avisynth.h>
#include <VapourSynth4.h>

struct DSVideoInfo
{
	DSFormat Format;
	int64_t FPSNum{ 1 };
	int64_t FPSDenom{ 1 };
	int Width{ 0 };
	int Height{ 0 };
	int Frames{ 0 };

	int Audio_SPS{ 0 };
	int Audio_SType{ 0 };
	int64_t Audio_NSamples{ 0 };
	int Audio_NChannels{ 0 };
	int Field{ 0 };

	DSVideoInfo() = default;

	DSVideoInfo(DSFormat format, int64_t fpsnum, int64_t fpsdenom, int width, int height, int frames)
		: Format(format)
		, FPSNum(fpsnum), FPSDenom(fpsdenom)
		, Width(width), Height(height)
		, Frames(frames)
	{
	}

	DSVideoInfo(const VSVideoInfo& vsvi)
		: Format(vsvi.format)
		, FPSNum(vsvi.fpsNum), FPSDenom(vsvi.fpsDen)
		, Width(vsvi.width), Height(vsvi.height)
		, Frames(vsvi.numFrames)
	{
	}

	DSVideoInfo(const VSVideoInfo* vsvi)
	{
		if (vsvi)
			*this = DSVideoInfo(*vsvi);
	}

	DSVideoInfo(const VideoInfo& avsvi)
		: Format(avsvi)
		, FPSNum(avsvi.fps_numerator), FPSDenom(avsvi.fps_denominator)
		, Width(avsvi.width), Height(avsvi.height)
		, Frames(avsvi.num_frames)
		, Audio_SPS(avsvi.audio_samples_per_second)
		, Audio_SType(avsvi.sample_type)
		, Audio_NSamples(avsvi.num_audio_samples)
		, Audio_NChannels(avsvi.nchannels)
		, Field(avsvi.image_type)
	{
	}

	VSVideoInfo ToVSVI(VSCore* core, const VSAPI* vsapi) const {
		VSVideoInfo vi{};
		vi.fpsNum = FPSNum;
		vi.fpsDen = FPSDenom;
		vi.width = Width;
		vi.height = Height;
		vi.numFrames = Frames;

		if (!Format.ToVSFormat(&vi.format, core, vsapi))
			throw "Unable to convert DSFormat to VSVideoFormat in ToVSVI";

		return vi;
	}

	VideoInfo ToAVSVI() const {
		VideoInfo vi{};
		vi.width = Width;
		vi.height = Height;
		vi.fps_numerator = static_cast<unsigned>(FPSNum);
		vi.fps_denominator = static_cast<unsigned>(FPSDenom);
		vi.num_frames = Frames;
		vi.pixel_type = Format.ToAVSFormat();
		vi.audio_samples_per_second = Audio_SPS;
		vi.sample_type = Audio_SType;
		vi.num_audio_samples = Audio_NSamples;
		vi.nchannels = Audio_NChannels;
		vi.image_type = Field;
		return vi;
	}
};
