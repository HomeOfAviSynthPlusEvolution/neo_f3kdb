#pragma once

#include "avisynth.h"

#include <stdio.h>

static void check_video_format(const char* name, const VideoInfo& vi, IScriptEnvironment* env)
{
    char name_buffer[256];
    if (!vi.IsYUV() || !vi.IsPlanar()) {
        sprintf_s(name_buffer, "%s: Only planar YUV clips are supported.", name);
        env->ThrowError(_strdup(name_buffer));
    }
    if (vi.IsFieldBased()) {
        sprintf_s(name_buffer, "%s: Field-based clip is not supported.", name);
        env->ThrowError(_strdup(name_buffer));
    }
}