/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - VapourSynth.
 *
 */

#pragma once

#include "ds_filter.hpp"
#include <cstdio>
#include <memory>
#include <vector>
#include <string>

namespace Plugin {
	extern const char* Identifier;
	extern const char* Namespace;
	extern const char* Description;
}

namespace VSInterface {

	inline thread_local VSFrameContext* tls_frameCtx{ nullptr };

	struct FrameContextGuard {
		VSFrameContext* prev_ctx;
		explicit FrameContextGuard(VSFrameContext* ctx) : prev_ctx(tls_frameCtx) {
			tls_frameCtx = ctx;
		}
		~FrameContextGuard() {
			tls_frameCtx = prev_ctx;
		}
		FrameContextGuard(const FrameContextGuard&) = delete;
		FrameContextGuard& operator=(const FrameContextGuard&) = delete;
	};

	struct VSInDelegator final : InDelegator {
		const VSMap* _in;
		const VSAPI* _vsapi;

		void* GetEnv() override { return nullptr; }
		bool IsAVS12() const override { return false; }

		void Read(const char* name, int& output) override {
			int err;
			int result = _vsapi->mapGetIntSaturated(_in, name, 0, &err);
			if (!err) output = result;
		}

		void Read(const char* name, int64_t& output) override {
			int err;
			int64_t result = _vsapi->mapGetInt(_in, name, 0, &err);
			if (!err) output = result;
		}

		void Read(const char* name, float& output) override {
			int err;
			float result = _vsapi->mapGetFloatSaturated(_in, name, 0, &err);
			if (!err) output = result;
		}

		void Read(const char* name, double& output) override {
			int err;
			double result = _vsapi->mapGetFloat(_in, name, 0, &err);
			if (!err) output = result;
		}

		void Read(const char* name, bool& output) override {
			int err;
			int64_t result = _vsapi->mapGetInt(_in, name, 0, &err);
			if (!err) output = (result != 0);
		}

		void Read(const char* name, std::string& output) override {
			int err;
			const char* result = _vsapi->mapGetData(_in, name, 0, &err);
			if (!err && result) {
				int size = _vsapi->mapGetDataSize(_in, name, 0, &err);
				if (!err && size >= 0)
					output.assign(result, static_cast<size_t>(size));
				else
					output = result;
			}
		}

		void Read(const char* name, std::vector<int>& output) override {
			int numElements = _vsapi->mapNumElements(_in, name);
			if (numElements <= 0) return;
			output.clear();
			output.reserve(numElements);
			for (int i = 0; i < numElements; i++) {
				int err;
				int val = _vsapi->mapGetIntSaturated(_in, name, i, &err);
				if (!err) output.push_back(val);
			}
		}

		void Read(const char* name, std::vector<int64_t>& output) override {
			int numElements = _vsapi->mapNumElements(_in, name);
			if (numElements <= 0) return;
			output.clear();
			output.reserve(numElements);
			for (int i = 0; i < numElements; i++) {
				int err;
				int64_t val = _vsapi->mapGetInt(_in, name, i, &err);
				if (!err) output.push_back(val);
			}
		}

		void Read(const char* name, std::vector<float>& output) override {
			int numElements = _vsapi->mapNumElements(_in, name);
			if (numElements <= 0) return;
			output.clear();
			output.reserve(numElements);
			for (int i = 0; i < numElements; i++) {
				int err;
				float val = _vsapi->mapGetFloatSaturated(_in, name, i, &err);
				if (!err) output.push_back(val);
			}
		}

		void Read(const char* name, std::vector<double>& output) override {
			int numElements = _vsapi->mapNumElements(_in, name);
			if (numElements <= 0) return;
			output.clear();
			output.reserve(numElements);
			for (int i = 0; i < numElements; i++) {
				int err;
				double val = _vsapi->mapGetFloat(_in, name, i, &err);
				if (!err) output.push_back(val);
			}
		}

		void Read(const char* name, std::vector<bool>& output) override {
			int numElements = _vsapi->mapNumElements(_in, name);
			if (numElements <= 0) return;
			output.clear();
			output.reserve(numElements);
			for (int i = 0; i < numElements; i++) {
				int err;
				int64_t val = _vsapi->mapGetInt(_in, name, i, &err);
				if (!err) output.push_back(val != 0);
			}
		}

		void Read(const char* name, void*& output) override {
			int err;
			VSNode* node = _vsapi->mapGetNode(_in, name, 0, &err);
			if (!err && node) {
				if (_vsapi->getNodeType(node) != mtVideo) {
					_vsapi->freeNode(node);
					throw "Clip parameter must be a video clip.";
				}
				output = reinterpret_cast<void*>(node);
			}
		}

		void Free(void*& clip) override {
			if (clip) {
				_vsapi->freeNode(reinterpret_cast<VSNode*>(clip));
				clip = nullptr;
			}
		}

		VSInDelegator(const VSMap* in, const VSAPI* vsapi)
			: _in(in), _vsapi(vsapi)
		{
		}
	};

	struct VSFetchFrameFunctor final : FetchFrameFunctor {
		VSNode* _vs_clip;
		VSCore* _core;
		const VSAPI* _vsapi;

		VSFetchFrameFunctor(VSNode* clip, VSCore* core, const VSAPI* vsapi)
			: _vs_clip(clip), _core(core), _vsapi(vsapi)
		{
		}

		DSFrame operator()(int n) override {
			if (!tls_frameCtx) {
				throw "DualSynth: fetch_frame cannot be called from unmanaged worker threads or outside GetFrame.";
			}

			const VSFrame* frame = _vsapi->getFrameFilter(n, _vs_clip, tls_frameCtx);
			if (!frame) {
				throw "DualSynth: Requested frame was not declared upfront in RequestReferenceFrames().";
			}

			return DSFrame(frame, _core, _vsapi);
		}

		~VSFetchFrameFunctor() override {
			if (_vs_clip) {
				_vsapi->freeNode(_vs_clip);
				_vs_clip = nullptr;
			}
		}
	};

	template<typename FilterType>
	void VS_CC Delete(void* instanceData, VSCore* core, const VSAPI* vsapi) {
		(void)core;
		(void)vsapi;
		auto filter = reinterpret_cast<FilterType*>(instanceData);
		if (filter) {
			if (filter->fetch_frame) {
				auto functor = reinterpret_cast<VSFetchFrameFunctor*>(filter->fetch_frame);
				delete functor;
				filter->fetch_frame = nullptr;
			}
			delete filter;
		}
	}

	template<typename FilterType>
	const VSFrame* VS_CC GetFrame(int n, int activationReason, void* instanceData,
		void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi) {
		(void)frameData;
		auto filter = reinterpret_cast<FilterType*>(instanceData);
		auto functor = reinterpret_cast<VSFetchFrameFunctor*>(filter->fetch_frame);

		try {
			if (activationReason == arInitial) {
				if (functor) {
					std::vector<int> ref_frames = filter->RequestReferenceFrames(n);
					if (ref_frames.empty()) {
						std::unordered_map<int, DSFrame> in_frames;
						FrameContextGuard ctx_guard(frameCtx);
						DSFrame out_frame = filter->GetFrame(n, in_frames);
						return out_frame.ToVSFrame();
					}
					for (auto&& i : ref_frames)
						vsapi->requestFrameFilter(i, functor->_vs_clip, frameCtx);
				}
				else {
					// Source filter: emit frame on arInitial
					std::unordered_map<int, DSFrame> in_frames;
					in_frames[n] = DSFrame(core, vsapi);

					FrameContextGuard ctx_guard(frameCtx);
					DSFrame out_frame = filter->GetFrame(n, in_frames);
					return out_frame.ToVSFrame();
				}
			}
			else if (activationReason == arAllFramesReady) {
				std::unordered_map<int, DSFrame> in_frames;

				if (functor) {
					std::vector<int> ref_frames = filter->RequestReferenceFrames(n);
					for (auto&& i : ref_frames)
						in_frames[i] = DSFrame(vsapi->getFrameFilter(i, functor->_vs_clip, frameCtx), core, vsapi);
				}
				else {
					in_frames[n] = DSFrame(core, vsapi);
				}

				FrameContextGuard ctx_guard(frameCtx);
				DSFrame out_frame = filter->GetFrame(n, in_frames);
				return out_frame.ToVSFrame();
			}
			else if (activationReason == arError) {
				return nullptr;
			}
		}
		catch (const char* err) {
			char msg_buff[512];
			std::snprintf(msg_buff, sizeof(msg_buff), "%s: %s", filter->VSName(), err);
			vsapi->setFilterError(msg_buff, frameCtx);
			return nullptr;
		}
		catch (const std::exception& e) {
			char msg_buff[512];
			std::snprintf(msg_buff, sizeof(msg_buff), "%s: %s", filter->VSName(), e.what());
			vsapi->setFilterError(msg_buff, frameCtx);
			return nullptr;
		}
		catch (...) {
			char msg_buff[512];
			std::snprintf(msg_buff, sizeof(msg_buff), "%s: Unknown exception in GetFrame", filter->VSName());
			vsapi->setFilterError(msg_buff, frameCtx);
			return nullptr;
		}

		return nullptr;
	}

	template<typename FilterType>
	void VS_CC Create(const VSMap* in, VSMap* out, void* userData, VSCore* core, const VSAPI* vsapi) {
		(void)userData;
		auto filter = std::make_unique<FilterType>();
		auto argument = VSInDelegator(in, vsapi);

		struct ClipNodesGuard {
			const VSAPI* vsapi;
			std::vector<VSNode*> nodes;
			size_t owned_from{ 0 };

			explicit ClipNodesGuard(const VSAPI* api) : vsapi(api) {}
			~ClipNodesGuard() {
				for (size_t i = owned_from; i < nodes.size(); ++i) {
					if (nodes[i])
						vsapi->freeNode(nodes[i]);
				}
			}
			ClipNodesGuard(const ClipNodesGuard&) = delete;
			ClipNodesGuard& operator=(const ClipNodesGuard&) = delete;
		} clip_guard(vsapi);

		std::unique_ptr<VSFetchFrameFunctor> functor;

		try {
			auto params = filter->Params();
			for (auto&& p : params) {
				if (p.Type == Clip && p.VSEnabled) {
					int err;
					VSNode* node = vsapi->mapGetNode(in, p.Name, 0, &err);
					if (!err && node) {
						if (vsapi->getNodeType(node) != mtVideo) {
							vsapi->freeNode(node);
							throw "Clip parameter must be a video clip.";
						}
						clip_guard.nodes.push_back(node);
					}
				}
			}

			DSVideoInfo input_vi;
			if (!clip_guard.nodes.empty()) {
				VSNode* primary_clip = clip_guard.nodes[0];
				const VSVideoInfo* vi_ptr = vsapi->getVideoInfo(primary_clip);
				if (!vi_ptr)
					throw "Unable to query VideoInfo from primary clip.";
				input_vi = DSVideoInfo(vi_ptr);

				functor = std::make_unique<VSFetchFrameFunctor>(primary_clip, core, vsapi);
				clip_guard.owned_from = 1;
			}

			filter->Initialize(&argument, input_vi, functor.get());

			auto output_vi = filter->GetOutputVI();
			VSVideoInfo vs_vi = output_vi.ToVSVI(core, vsapi);

			std::vector<VSFilterDependency> dependencies;
			dependencies.reserve(clip_guard.nodes.size());
			for (auto node : clip_guard.nodes) {
				dependencies.push_back({ node, filter->GetVSRequestPattern() });
			}

			VSNode* node = vsapi->createVideoFilter2(
				filter->VSName(),
				&vs_vi,
				GetFrame<FilterType>,
				Delete<FilterType>,
				filter->VSMode(),
				dependencies.data(),
				static_cast<int>(dependencies.size()),
				filter.get(),
				core
			);

			if (!node) {
				vsapi->mapSetError(out, "Failed to create filter node");
				return;
			}

			// Ownership transferred to VapourSynth core
			functor.release();
			filter.release();

			// clip_guard destructor runs on scope exit, releasing clip_guard.nodes[1...N]
			vsapi->mapConsumeNode(out, "clip", node, maAppend);
		}
		catch (const char* err) {
			char msg_buff[512];
			std::snprintf(msg_buff, sizeof(msg_buff), "%s: %s", filter->VSName(), err);
			vsapi->mapSetError(out, msg_buff);
		}
		catch (const std::exception& e) {
			char msg_buff[512];
			std::snprintf(msg_buff, sizeof(msg_buff), "%s: %s", filter->VSName(), e.what());
			vsapi->mapSetError(out, msg_buff);
		}
		catch (...) {
			char msg_buff[512];
			std::snprintf(msg_buff, sizeof(msg_buff), "%s: Unknown exception during filter creation", filter->VSName());
			vsapi->mapSetError(out, msg_buff);
		}
	}

	template<typename FilterType>
	void RegisterFilter(VSPlugin* plugin, const VSPLUGINAPI* vspapi) {
		FilterType filter;
		int result = vspapi->registerFunction(
			filter.VSName(),
			filter.VSParams().c_str(),
			filter.VSReturnType().c_str(),
			Create<FilterType>,
			nullptr,
			plugin
		);
		(void)result;
	}

} // namespace VSInterface

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin* plugin, const VSPLUGINAPI* vspapi) {
	vspapi->configPlugin(
		Plugin::Identifier,
		Plugin::Namespace,
		Plugin::Description,
		VS_MAKE_VERSION(1, 0),
		VAPOURSYNTH_API_VERSION,
		0,
		plugin
	);

	auto filters = RegisterVSFilters();
	for (auto&& RegisterFilter : filters) {
		RegisterFilter(plugin, vspapi);
	}
}
