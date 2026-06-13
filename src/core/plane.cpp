#include "core/plane.hpp"

#include <cstring>

void destroy_context(process_plane_context* context) {
  if (context->destroy != nullptr) {
    context->destroy(context->data);
  }
  std::memset(context, 0, sizeof(process_plane_context));
}

void init_context(process_plane_context* context) {
  std::memset(context, 0, sizeof(process_plane_context));
}
