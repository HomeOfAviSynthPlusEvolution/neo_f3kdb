/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - Common header+.
 *
 */

#pragma once

#include <avisynth.h>
#include <VapourSynth4.h>
#include <cstring>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <mutex>
#include <functional>

#include "ds_format.hpp"
#include "ds_videoinfo.hpp"
#include "ds_frame.hpp"

inline std::mutex& GetFFTWMutex() {
    static std::mutex m;
    return m;
}

class GlobalLockGuard
{
    IScriptEnvironment* env;
    const char* name;
    bool acquired;
    bool is_legacy;

public:
    GlobalLockGuard(IScriptEnvironment* _env, const char* _name, bool use_avs_lock)
        : env(_env), name(_name), acquired(false), is_legacy(false)
    {
        if (!name) return;

        if (env && use_avs_lock) {
            acquired = env->AcquireGlobalLock(name);
            if (acquired)
                return;

            throw std::runtime_error("Failed to acquire AviSynth global lock: " + std::string(name));
        }

        if (std::strcmp(name, "fftw") == 0) {
            GetFFTWMutex().lock();
            acquired = true;
            is_legacy = true;
        }
        else {
            throw std::runtime_error("No global lock mechanism available for: " + std::string(name));
        }
    }

    ~GlobalLockGuard() {
        if (acquired) {
            if (is_legacy)
                GetFFTWMutex().unlock();
            else if (env)
                env->ReleaseGlobalLock(name);
        }
    }

    GlobalLockGuard(const GlobalLockGuard&) = delete;
    GlobalLockGuard& operator=(const GlobalLockGuard&) = delete;
};

typedef void (*register_vsfilter_proc)(VSPlugin* plugin, const VSPLUGINAPI* vspapi);

typedef void (*register_avsfilter_proc)(IScriptEnvironment* env);

std::vector<register_vsfilter_proc> RegisterVSFilters();
std::vector<register_avsfilter_proc> RegisterAVSFilters();

enum ParamType
{
	Clip,
	Integer,
	Float,
	Boolean,
	String
};

struct Param
{
    const char* Name{ nullptr };
    ParamType Type{ Integer };
    bool IsArray{ false };
    bool AVSEnabled{ true };
    bool VSEnabled{ true };
    bool IsOptional{ true };
};

struct InDelegator
{
    virtual void* GetEnv() { return nullptr; }
    virtual bool IsAVS12() const { return false; }

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

	virtual ~InDelegator() = default;
};

struct FetchFrameFunctor
{
	virtual DSFrame operator()(int n) = 0;
	virtual ~FetchFrameFunctor() = default;
};
