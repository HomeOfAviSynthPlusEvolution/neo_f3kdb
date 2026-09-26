/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - DSFormat.
 *
 */

#pragma once

struct DSFormat
{
	bool IsFamilyYUV{ true };
	bool IsFamilyRGB{ false };
	bool IsFamilyGray{ false };
	bool IsInteger{ true };
	bool IsFloat{ false };
	int SSW{ 0 };
	int SSH{ 0 };
	int BitsPerSample{ 8 };
	int BytesPerSample{ 1 };
	int Planes{ 3 };

	DSFormat() = default;

	DSFormat(const VSVideoFormat& format)
	{
		Planes = format.numPlanes;
		IsFamilyYUV = (format.colorFamily == cfYUV);
		IsFamilyRGB = (format.colorFamily == cfRGB);
		IsFamilyGray = (format.colorFamily == cfGray);
		SSW = format.subSamplingW;
		SSH = format.subSamplingH;
		BitsPerSample = format.bitsPerSample;
		BytesPerSample = format.bytesPerSample;
		IsInteger = (format.sampleType == stInteger);
		IsFloat = (format.sampleType == stFloat);
	}

	bool ToVSFormat(VSVideoFormat* outFormat, VSCore* core, const VSAPI* vsapi) const
	{
		if (!outFormat || !vsapi || Planes > 3)
			return false;

		int colorFamily = cfUndefined;
		if (IsFamilyYUV)
			colorFamily = cfYUV;
		else if (IsFamilyRGB)
			colorFamily = cfRGB;
		else if (IsFamilyGray)
			colorFamily = cfGray;

		int sampleType = IsFloat ? stFloat : stInteger;
		// VapourSynth requires subSamplingW and subSamplingH to be 0 for RGB and Gray formats
		int ssw = (colorFamily == cfYUV) ? SSW : 0;
		int ssh = (colorFamily == cfYUV) ? SSH : 0;

		return vsapi->queryVideoFormat(outFormat, colorFamily, sampleType, BitsPerSample, ssw, ssh, core) != 0;
	}

	DSFormat(const VideoInfo& vi)
	{
		if (!vi.IsPlanar())
			throw "DualSynth only supports planar formats.";

		IsFamilyGray = vi.IsY();
		IsFamilyYUV = (vi.IsYUV() || vi.IsYUVA()) && !IsFamilyGray;
		IsFamilyRGB = vi.IsRGB() || vi.IsPlanarRGBA();

		IsFloat = (vi.ComponentSize() == 4);
		IsInteger = !IsFloat;
		Planes = vi.NumComponents();
		BitsPerSample = vi.BitsPerComponent();
		BytesPerSample = vi.ComponentSize();

		if (IsFamilyYUV && Planes > 1) {
			SSW = vi.GetPlaneWidthSubsampling(PLANAR_U);
			SSH = vi.GetPlaneHeightSubsampling(PLANAR_U);
		}
	}

	int ToAVSFormat() const
	{
		int pixel_format = 0;
		if (IsFamilyGray) {
			pixel_format = VideoInfo::CS_GENERIC_Y;
		}
		else if (IsFamilyYUV) {
			if (Planes == 1) {
				pixel_format = VideoInfo::CS_GENERIC_Y;
			}
			else {
				pixel_format = VideoInfo::CS_PLANAR | (Planes == 4 ? VideoInfo::CS_YUVA : VideoInfo::CS_YUV) | VideoInfo::CS_VPlaneFirst;
				switch (SSW) {
				case 0: pixel_format |= VideoInfo::CS_Sub_Width_1; break;
				case 1: pixel_format |= VideoInfo::CS_Sub_Width_2; break;
				case 2: pixel_format |= VideoInfo::CS_Sub_Width_4; break;
				}
				switch (SSH) {
				case 0: pixel_format |= VideoInfo::CS_Sub_Height_1; break;
				case 1: pixel_format |= VideoInfo::CS_Sub_Height_2; break;
				case 2: pixel_format |= VideoInfo::CS_Sub_Height_4; break;
				}
			}
		}
		else if (IsFamilyRGB) {
			pixel_format = VideoInfo::CS_PLANAR | VideoInfo::CS_BGR | (Planes == 4 ? VideoInfo::CS_RGBA_TYPE : VideoInfo::CS_RGB_TYPE);
		}

		switch (BitsPerSample) {
		case 8: pixel_format |= VideoInfo::CS_Sample_Bits_8; break;
		case 10: pixel_format |= VideoInfo::CS_Sample_Bits_10; break;
		case 12: pixel_format |= VideoInfo::CS_Sample_Bits_12; break;
		case 14: pixel_format |= VideoInfo::CS_Sample_Bits_14; break;
		case 16: pixel_format |= VideoInfo::CS_Sample_Bits_16; break;
		case 32: pixel_format |= VideoInfo::CS_Sample_Bits_32; break;
		}
		return pixel_format;
	}
};
