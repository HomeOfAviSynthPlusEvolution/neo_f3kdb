#pragma once

#include "core/constants.hpp"
#include "core/kernel_types.hpp"

#define DEFAULT_RANDOM_PARAM 1.0

int random(RANDOM_ALGORITHM algo, int& seed, int range, double param);
