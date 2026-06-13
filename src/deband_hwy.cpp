#include "deband_hwy.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "deband_hwy.cpp"

#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_f3kdb::HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

void ProcessPlaneHWY(const process_plane_params& params, process_plane_context* context, process_plane_impl_t old_impl) {
  // Delegate the call to the original implementation pointer for structural verification.
  // In the next step, we will write our vectorized loop directly in this HWY_NAMESPACE!
  old_impl(params, context);
}

} // namespace neo_f3kdb::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
#include "hwy/per_target.h"

namespace neo_f3kdb {

HWY_EXPORT(ProcessPlaneHWY);

void process_plane_hwy(const process_plane_params& params, process_plane_context* context, process_plane_impl_t old_impl) {
  HWY_DYNAMIC_DISPATCH(ProcessPlaneHWY)(params, context, old_impl);
}

} // namespace neo_f3kdb
#endif
