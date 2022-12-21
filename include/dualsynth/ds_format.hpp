/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - DSFormat.
 *
 */

#pragma once

struct DSFormat
{
  bool IsFamilyYUV {true}, IsFamilyRGB {false}, IsFamilyYCC {false};
  bool IsInteger {true}, IsFloat {false};
  int SSW {0}, SSH {0};
  int BitsPerSample {8}, BytesPerSample {1};
  int Planes {3};
  DSFormat() {}
  DSFormat(const VSVideoFormat format)
  {
    Planes = format.numPlanes;
    IsFamilyYUV = format.colorFamily == cfYUV || format.colorFamily == cfGray;
    IsFamilyRGB = format.colorFamily == cfRGB;
    // IsFamilyYCC = format.colorFamily == cfYCoCg;
    SSW = format.subSamplingW;
    SSH = format.subSamplingH;
    BitsPerSample = format.bitsPerSample;
    BytesPerSample = format.bytesPerSample;
    IsInteger = format.sampleType == stInteger;
    IsFloat = format.sampleType == stFloat;
  }

  const VSVideoFormat ToVSFormat() const
  {
    VSColorFamily family = cfYUV;
    if (IsFamilyYUV)
      family = Planes == 1 ? cfGray : cfYUV;
    else if (IsFamilyRGB)
      family = cfRGB;
    // else if (IsFamilyYCC)
    //   family = cfYCoCg;
    
    VSVideoFormat format;
    format.colorFamily = family;
    format.sampleType = IsInteger ? stInteger : stFloat;
    format.bitsPerSample = BitsPerSample;
    format.bytesPerSample = BytesPerSample;
    format.subSamplingW = SSW;
    format.subSamplingH = SSH;
    format.numPlanes = Planes;
    return format;
  }

  DSFormat(int format)
  {
    const int componentBitSizes[8] = {8,16,32,0,0,10,12,14};
    if (format == VideoInfo::CS_I420)
      format = VideoInfo::CS_YV12;

    auto PYUV = VideoInfo::CS_PLANAR | VideoInfo::CS_YUV;
    IsFamilyYUV = (format & PYUV) == PYUV;
    auto PRGB = VideoInfo::CS_PLANAR | VideoInfo::CS_BGR;
    IsFamilyRGB = (format & PRGB) == PRGB;
    IsFamilyYCC = false;
    BitsPerSample = componentBitSizes[(format >> VideoInfo::CS_Shift_Sample_Bits) & 7];
    BytesPerSample = BitsPerSample == 8 ? 1 : BitsPerSample == 32 ? 4 : 2;
    IsInteger = BitsPerSample < 32;
    IsFloat = BitsPerSample == 32;
    if (IsFamilyYUV && (format & VideoInfo::CS_GENERIC_Y) == VideoInfo::CS_GENERIC_Y)
      Planes = 1;
    else if (IsFamilyYUV && (format & VideoInfo::CS_YUVA) == VideoInfo::CS_YUVA)
      Planes = 4;
    else if (IsFamilyRGB && (format & VideoInfo::CS_RGBA_TYPE) == VideoInfo::CS_RGBA_TYPE)
      Planes = 4;

    if (IsFamilyYUV && Planes > 1) {
      SSW = ((format >> VideoInfo::CS_Shift_Sub_Width) + 1) & 3;
      SSH = ((format >> VideoInfo::CS_Shift_Sub_Height) + 1) & 3;
    }
  }

  int ToAVSFormat() const
  {
    int pixel_format = VideoInfo::CS_PLANAR | (Planes == 3 ? VideoInfo::CS_YUV : VideoInfo::CS_YUVA) | VideoInfo::CS_VPlaneFirst;
    if (IsFamilyYUV) {
      pixel_format = VideoInfo::CS_PLANAR | (Planes == 3 ? VideoInfo::CS_YUV : VideoInfo::CS_YUVA) | VideoInfo::CS_VPlaneFirst;

      switch(SSW) {
        case 0: pixel_format |= VideoInfo::CS_Sub_Width_1; break;
        case 1: pixel_format |= VideoInfo::CS_Sub_Width_2; break;
        case 2: pixel_format |= VideoInfo::CS_Sub_Width_4; break;
      }

      switch(SSH) {
        case 0: pixel_format |= VideoInfo::CS_Sub_Height_1; break;
        case 1: pixel_format |= VideoInfo::CS_Sub_Height_2; break;
        case 2: pixel_format |= VideoInfo::CS_Sub_Height_4; break;
      }

      if (Planes == 1)
        pixel_format = VideoInfo::CS_GENERIC_Y;
    }
    else if (IsFamilyRGB || IsFamilyYCC)
      pixel_format = VideoInfo::CS_PLANAR | VideoInfo::CS_BGR | (Planes == 3 ? VideoInfo::CS_RGB_TYPE : VideoInfo::CS_RGBA_TYPE);

    switch(BitsPerSample) {
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
