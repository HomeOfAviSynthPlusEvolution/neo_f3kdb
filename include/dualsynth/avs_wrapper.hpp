/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - AviSynth+.
 *
 */

#pragma once

#include "ds_filter.hpp"
#include <avisynth.h>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <mutex>

namespace Plugin {
	extern const char* Description;
}

namespace AVSInterface
{
	struct AVSInDelegator final : InDelegator {
		const AVSValue _args;
		IScriptEnvironment* _env;
		bool _is_v12;
		std::unordered_map<std::string, int> _params_index_map;

		int NameToIndex(const char* name) {
			std::string name_string(name ? name : "");
			auto it = _params_index_map.find(name_string);
			if (it == _params_index_map.end())
				throw "Unknown parameter during NameToIndex";
			return it->second;
		}

		void Read(const char* name, int& output) override {
			auto arg = _args[NameToIndex(name)];
			if (arg.Defined())
				output = arg.AsInt(output);
		}

		void Read(const char* name, int64_t& output) override {
			auto arg = _args[NameToIndex(name)];
			if (arg.Defined())
				output = arg.AsInt(output);
		}

		void Read(const char* name, float& output) override {
			auto arg = _args[NameToIndex(name)];
			if (arg.Defined())
				output = static_cast<float>(arg.AsFloat(output));
		}

		void Read(const char* name, double& output) override {
			auto arg = _args[NameToIndex(name)];
			if (arg.Defined())
				output = arg.AsFloat(output);
		}

		void Read(const char* name, bool& output) override {
			auto arg = _args[NameToIndex(name)];
			if (arg.Defined())
				output = arg.AsBool(output);
		}

		void Read(const char* name, std::string& output) override {
			auto arg = _args[NameToIndex(name)];
			if (arg.Defined()) {
				const char* result = arg.AsString(nullptr);
				if (result)
					output = result;
			}
		}

		void Read(const char* name, void*& output) override {
			auto arg = _args[NameToIndex(name)];
			if (arg.Defined() && arg.IsClip()) {
				PClip* clip = new PClip(arg.AsClip());
				output = reinterpret_cast<void*>(clip);
			}
		}

		void Read(const char* name, std::vector<int>& output) override {
			auto arg = _args[NameToIndex(name)];
			if (!arg.Defined())
				return;
			if (!arg.IsArray())
				throw "Argument is not array";
			auto size = arg.ArraySize();
			output.clear();
			output.reserve(size);
			for (int i = 0; i < size; i++)
				output.push_back(arg[i].AsInt());
		}

		void Read(const char* name, std::vector<int64_t>& output) override {
			auto arg = _args[NameToIndex(name)];
			if (!arg.Defined())
				return;
			if (!arg.IsArray())
				throw "Argument is not array";
			auto size = arg.ArraySize();
			output.clear();
			output.reserve(size);
			for (int i = 0; i < size; i++)
				output.push_back(arg[i].AsInt());
		}

		void Read(const char* name, std::vector<float>& output) override {
			auto arg = _args[NameToIndex(name)];
			if (!arg.Defined())
				return;
			if (!arg.IsArray())
				throw "Argument is not array";
			auto size = arg.ArraySize();
			output.clear();
			output.reserve(size);
			for (int i = 0; i < size; i++)
				output.push_back(static_cast<float>(arg[i].AsFloat()));
		}

		void Read(const char* name, std::vector<double>& output) override {
			auto arg = _args[NameToIndex(name)];
			if (!arg.Defined())
				return;
			if (!arg.IsArray())
				throw "Argument is not array";
			auto size = arg.ArraySize();
			output.clear();
			output.reserve(size);
			for (int i = 0; i < size; i++)
				output.push_back(arg[i].AsFloat());
		}

		void Read(const char* name, std::vector<bool>& output) override {
			auto arg = _args[NameToIndex(name)];
			if (!arg.Defined())
				return;
			if (!arg.IsArray())
				throw "Argument is not array";
			auto size = arg.ArraySize();
			output.clear();
			output.reserve(size);
			for (int i = 0; i < size; i++)
				output.push_back(arg[i].AsBool());
		}

		void Free(void*& clip) override {
			if (clip) {
				PClip* c = reinterpret_cast<PClip*>(clip);
				delete c;
				clip = nullptr;
			}
		}

		void* GetEnv() override { return _env; }
		bool IsAVS12() const override { return _is_v12; }

		AVSInDelegator(const AVSValue args, std::vector<Param> params, IScriptEnvironment* env)
			: _args(args), _env(env), _is_v12(false)
		{
			if (_env) {
				try {
					_env->CheckVersion(12);
					_is_v12 = true;
				}
				catch (...) {
					_is_v12 = false;
				}
			}

			int idx = 0;
			for (auto&& param : params)
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

		AVSFetchFrameFunctor(PClip clip, VideoInfo vi, IScriptEnvironment* env)
			: _clip(clip), _vi(vi), _env(env)
		{
		}

		DSFrame operator()(int n) override {
			std::lock_guard<std::mutex> guard(fetch_frame_mutex);
			auto frame = _clip->GetFrame(n, _env);
			return DSFrame(frame, _vi, _env);
		}

		~AVSFetchFrameFunctor() override = default;
	};

	template<typename FilterType>
	struct AVSWrapper : IClip
	{
		AVSValue _args;
		IScriptEnvironment* _env;
		FilterType data;
		PClip clip;
		VideoInfo vi;
		AVSFetchFrameFunctor* functor{ nullptr };

		AVSWrapper(AVSValue args, IScriptEnvironment* env)
			: _args(args), _env(env)
		{
		}

		void Initialize()
		{
			auto input_vi = DSVideoInfo();
			if (_args[0].IsClip()) {
				clip = _args[0].AsClip();
				input_vi = DSVideoInfo(clip->GetVideoInfo());
				functor = new AVSFetchFrameFunctor(clip, clip->GetVideoInfo(), _env);
			}
			auto argument = AVSInDelegator(_args, data.AVSOrderedParams(), _env);
			data.Initialize(&argument, input_vi, functor);
		}

		PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override {
			try {
				std::unordered_map<int, DSFrame> in_frames;
				if (functor) {
					std::vector<int> requests = data.RequestReferenceFrames(n);
					auto input_vi = clip->GetVideoInfo();
					for (auto&& i : requests) {
						auto frame = clip->GetFrame(i, env);
						in_frames[i] = DSFrame(frame, input_vi, env);
					}
				}
				else {
					in_frames[n] = DSFrame(env);
				}

				return data.GetFrame(n, in_frames).ToAVSFrame();
			}
			catch (const AvisynthError&) {
				throw;
			}
			catch (const std::exception& e) {
				env->ThrowError("%s: %s", data.AVSName(), e.what());
				return nullptr;
			}
			catch (const char* err) {
				env->ThrowError("%s: %s", data.AVSName(), err ? err : "Unknown error");
				return nullptr;
			}
			catch (...) {
				env->ThrowError("%s: Unknown exception in GetFrame", data.AVSName());
				return nullptr;
			}
		}

		const VideoInfo& __stdcall GetVideoInfo() override {
			auto output_vi = data.GetOutputVI();
			vi = output_vi.ToAVSVI();
			return vi;
		}

		void __stdcall GetAudio(void* buf, int64_t start, int64_t count, IScriptEnvironment* env) override {
			if (clip) clip->GetAudio(buf, start, count, env);
		}

		bool __stdcall GetParity(int n) override {
			return clip ? clip->GetParity(n) : false;
		}

		int __stdcall SetCacheHints(int cachehints, int frame_range) override {
			return data.SetCacheHints(cachehints, frame_range);
		}

		~AVSWrapper() {
			delete functor;
		}
	};

	template<typename FilterType>
	AVSValue __cdecl Create(AVSValue args, void* user_data, IScriptEnvironment* env)
	{
		(void)user_data;
		auto filter = std::make_unique<AVSWrapper<FilterType>>(args, env);
		std::string filter_name = filter->data.AVSName();

		try {
			filter->Initialize();
		}
		catch (const char* err) {
			filter.reset();
			env->ThrowError("%s: %s", filter_name.c_str(), err ? err : "Unknown error");
		}
		catch (const AvisynthError& e) {
			std::string msg = e.msg ? e.msg : "Unknown AvisynthError";
			filter.reset();
			env->ThrowError("%s: %s", filter_name.c_str(), msg.c_str());
		}
		catch (const std::exception& e) {
			std::string msg = e.what();
			filter.reset();
			env->ThrowError("%s: %s", filter_name.c_str(), msg.c_str());
		}
		catch (...) {
			filter.reset();
			env->ThrowError("%s: Unknown exception during initialization", filter_name.c_str());
		}

		return filter.release();
	}

	template<typename FilterType>
	void RegisterFilter(IScriptEnvironment* env) {
		FilterType filter;
		env->AddFunction(filter.AVSName(), filter.AVSParams().c_str(), Create<FilterType>, nullptr);
	}
}

const AVS_Linkage* AVS_linkage = nullptr;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, AVS_Linkage* linkage)
{
	AVS_linkage = linkage;
	auto filters = RegisterAVSFilters();
	for (auto&& RegisterFilter : filters) {
		RegisterFilter(env);
	}
	return Plugin::Description;
}
