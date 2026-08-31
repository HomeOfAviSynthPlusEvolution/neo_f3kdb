#pragma once

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <dualsynth/span2d.hpp>
#include <stdexcept>
#include <utility>

#if defined(_WIN32)
#include <malloc.h>
#endif

namespace neo_f3kdb {

template <class T, std::size_t Alignment>
class AlignedBuffer {
public:
  AlignedBuffer() = default;

  explicit AlignedBuffer(std::size_t count) {
    reset(count);
  }

  ~AlignedBuffer() {
    release();
  }

  AlignedBuffer(const AlignedBuffer&) = delete;
  AlignedBuffer& operator=(const AlignedBuffer&) = delete;

  AlignedBuffer(AlignedBuffer&& other) noexcept
    : data_(std::exchange(other.data_, nullptr)),
      count_(std::exchange(other.count_, 0)) {}

  AlignedBuffer& operator=(AlignedBuffer&& other) noexcept {
    if (this != &other) {
      release();
      data_ = std::exchange(other.data_, nullptr);
      count_ = std::exchange(other.count_, 0);
    }
    return *this;
  }

  void reset(std::size_t count) {
    release();
    count_ = count;
    if (count_ == 0) {
      return;
    }

    const std::size_t bytes = count_ * sizeof(T);
#if defined(_WIN32)
    data_ = static_cast<T*>(_aligned_malloc(bytes, Alignment));
#else
    void* ptr = nullptr;
    constexpr std::size_t min_alignment = sizeof(void*);
    constexpr std::size_t alloc_alignment = Alignment < min_alignment ? min_alignment : Alignment;
    if (posix_memalign(&ptr, alloc_alignment, bytes) != 0) {
      ptr = nullptr;
    }
    data_ = static_cast<T*>(ptr);
#endif
    if (data_ == nullptr) {
      throw std::runtime_error("unable to allocate memory");
    }
    std::memset(data_, 0, bytes);
  }

  [[nodiscard]] T* data() noexcept {
    return data_;
  }

  [[nodiscard]] const T* data() const noexcept {
    return data_;
  }

  [[nodiscard]] span2d::Span<T> span() noexcept {
    return {data_, count_};
  }

  [[nodiscard]] span2d::ReadOnlySpan<T> span() const noexcept {
    return {data_, count_};
  }

  [[nodiscard]] span2d::Span<T> subspan(std::size_t offset, std::size_t count) noexcept {
    return span().subspan(offset, count);
  }

  [[nodiscard]] span2d::ReadOnlySpan<T> subspan(std::size_t offset, std::size_t count) const noexcept {
    return span().subspan(offset, count);
  }

  [[nodiscard]] std::size_t size() const noexcept {
    return count_;
  }

  [[nodiscard]] bool empty() const noexcept {
    return count_ == 0;
  }

  T& operator[](std::size_t index) noexcept {
    return data_[index];
  }

  const T& operator[](std::size_t index) const noexcept {
    return data_[index];
  }

private:
  void release() noexcept {
    if (data_ == nullptr) {
      return;
    }
#if defined(_WIN32)
    _aligned_free(data_);
#else
    std::free(data_);
#endif
    data_ = nullptr;
    count_ = 0;
  }

  T* data_ = nullptr;
  std::size_t count_ = 0;
};

} // namespace neo_f3kdb
