/*
 * Copyright 2020 Xinyue Lu
 *
 * DualSynth bridge - plugin.
 *
 */

#pragma once

#include "version.hpp"
#include "f3kdb.hpp"

namespace Plugin {
  const char* Identifier = "in.7086.neo_f3kdb";
  const char* Namespace = "neo_f3kdb";
  const char* Description = "Neo F3KDB Deband Filter " PLUGIN_VERSION;
  const int Version = 8;
}

std::vector<register_vsfilter_proc> RegisterVSFilters()
{
  return std::vector<register_vsfilter_proc> { VSInterface::RegisterFilter<F3KDB> };
}

std::vector<register_avsfilter_proc> RegisterAVSFilters()
{
  return std::vector<register_avsfilter_proc> { AVSInterface::RegisterFilter<F3KDB> };
}
