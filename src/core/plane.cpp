#include "core/plane.hpp"

void destroy_context(process_plane_context* context) {
  if (context == nullptr) {
    return;
  }
  std::lock_guard<std::mutex> lock(context->mutex);
  void* ptr = context->data.load(std::memory_order_relaxed);
  if (context->destroy != nullptr && ptr != nullptr) {
    context->destroy(ptr);
  }
  context->data.store(nullptr, std::memory_order_relaxed);
  context->destroy = nullptr;
}

void init_context(process_plane_context* context) {
  if (context == nullptr) {
    return;
  }
  std::lock_guard<std::mutex> lock(context->mutex);
  void* ptr = context->data.load(std::memory_order_relaxed);
  if (context->destroy != nullptr && ptr != nullptr) {
    context->destroy(ptr);
  }
  context->data.store(nullptr, std::memory_order_relaxed);
  context->destroy = nullptr;
}
