/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth wrapper - Filter parent class.
 *
 */

#pragma once

#include "ds_common.hpp"

struct Filter
{
	DSVideoInfo in_vi;
	FetchFrameFunctor* fetch_frame{ nullptr };

	virtual const char* VSName() const { return "FilterFoo"; }
	virtual const char* AVSName() const { return "FilterFoo"; }
	virtual MtMode AVSMode() const { return MT_SERIALIZED; }

	virtual VSFilterMode VSMode() const { return fmParallelRequests; }

	virtual VSRequestPattern GetVSRequestPattern() const { return rpGeneral; }

	virtual std::vector<Param> Params() const = 0;

	virtual std::vector<Param> AVSOrderedParams() const
	{
		auto params = this->Params();
		std::vector<Param> ordered;
		ordered.reserve(params.size());
		for (const auto& p : params) {
			if (p.AVSEnabled && !p.IsArray)
				ordered.push_back(p);
		}
		for (const auto& p : params) {
			if (p.AVSEnabled && p.IsArray)
				ordered.push_back(p);
		}
		return ordered;
	}

	virtual std::string VSParams() const
	{
		std::stringstream ss;
		auto params = this->Params();
		for (auto&& p : params)
		{
			if (!p.VSEnabled) continue;

			std::string type_name;
			switch (p.Type) {
			case Clip:    type_name = "vnode"; break;
			case Integer: type_name = "int";   break;
			case Float:   type_name = "float"; break;
			case Boolean: type_name = "int";   break;
			case String:  type_name = "data";  break;
			default:      type_name = "int";   break;
			}

			ss << p.Name << ':' << type_name;
			if (p.IsArray)
				ss << "[]";
			if (p.IsOptional)
				ss << ":opt";
			ss << ';';
		}
		return ss.str();
	}

	virtual std::string VSReturnType() const
	{
		return "clip:vnode;";
	}

	virtual std::string AVSParams() const
	{
		std::stringstream ss;
		auto params = this->AVSOrderedParams();
		for (auto&& p : params)
		{
			char type_name;
			switch (p.Type) {
			case Clip:    type_name = 'c'; break;
			case Integer: type_name = 'i'; break;
			case Float:   type_name = 'f'; break;
			case Boolean: type_name = 'b'; break;
			case String:  type_name = 's'; break;
			default:      type_name = 'i'; break;
			}

			if (p.IsArray) {
				ss << '[' << p.Name << ']' << type_name << (p.IsOptional ? '*' : '+');
			}
			else {
				if (p.IsOptional) ss << '[' << p.Name << ']';
				ss << type_name;
			}
		}
		return ss.str();
	}

	virtual void Initialize(InDelegator* in, DSVideoInfo in_vi, FetchFrameFunctor* fetch_frame)
	{
		(void)in;
		this->in_vi = in_vi;
		this->fetch_frame = fetch_frame;
	}

	virtual std::vector<int> RequestReferenceFrames(int n) const
	{
		return std::vector<int>{n};
	}

	virtual DSFrame GetFrame(int n, const std::unordered_map<int, DSFrame>& in_frames)
	{
		(void)n;
		return in_frames.empty() ? DSFrame() : in_frames.begin()->second;
	}

	virtual DSVideoInfo GetOutputVI()
	{
		return in_vi;
	}

	virtual int SetCacheHints(int cachehints, int frame_range)
	{
		(void)frame_range;
		return cachehints == CACHE_GET_MTMODE ? AVSMode() : 0;
	}

	virtual ~Filter() = default;
};
