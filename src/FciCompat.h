/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 */

/* Compatibility helpers used to keep the active-space FCI kernel independent
 * of external Boost libraries.
 */
#ifndef FCICOMPAT_H
#define FCICOMPAT_H

#include <cstdint>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace boost {
using uint64_t = std::uint64_t;
using uint32_t = std::uint32_t;
using uint16_t = std::uint16_t;
using uint8_t = std::uint8_t;
using int64_t = std::int64_t;
using int32_t = std::int32_t;
using int16_t = std::int16_t;
using int8_t = std::int8_t;

class format {
public:
  explicit format(std::string fmt) : fmt_(std::move(fmt)) {}

  template <typename T> format &operator%(const T &value) {
    std::ostringstream out;
    out << value;
    args_.push_back(out.str());
    return *this;
  }

  std::string str() const {
    std::ostringstream out;
    std::size_t arg = 0;
    for (std::size_t i = 0; i < fmt_.size(); ++i) {
      if (fmt_[i] != '%' || i + 1 == fmt_.size()) {
        out << fmt_[i];
        continue;
      }
      if (fmt_[i + 1] == '%') {
        out << '%';
        ++i;
        continue;
      }
      std::size_t j = i + 1;
      while (j < fmt_.size() &&
             std::string("diuoxXfFeEgGaAcsp").find(fmt_[j]) ==
                 std::string::npos) {
        ++j;
      }
      if (arg < args_.size()) {
        out << args_[arg++];
      }
      i = (j < fmt_.size()) ? j : fmt_.size() - 1;
    }
    return out.str();
  }

private:
  std::string fmt_;
  std::vector<std::string> args_;
};

inline std::string str(const format &fmt) { return fmt.str(); }

inline std::ostream &operator<<(std::ostream &out, const format &fmt) {
  out << fmt.str();
  return out;
}

template <typename T> class intrusive_ptr {
public:
  intrusive_ptr() : ptr_(nullptr) {}
  intrusive_ptr(T *ptr, bool add_ref = true) : ptr_(ptr) {
    if (ptr_ != nullptr && add_ref)
      intrusive_ptr_add_ref(ptr_);
  }
  intrusive_ptr(const intrusive_ptr &other) : ptr_(other.ptr_) {
    if (ptr_ != nullptr)
      intrusive_ptr_add_ref(ptr_);
  }
  intrusive_ptr(intrusive_ptr &&other) noexcept : ptr_(other.ptr_) {
    other.ptr_ = nullptr;
  }
  ~intrusive_ptr() {
    if (ptr_ != nullptr)
      intrusive_ptr_release(ptr_);
  }

  intrusive_ptr &operator=(const intrusive_ptr &other) {
    reset(other.ptr_);
    return *this;
  }
  intrusive_ptr &operator=(intrusive_ptr &&other) noexcept {
    if (this != &other) {
      if (ptr_ != nullptr)
        intrusive_ptr_release(ptr_);
      ptr_ = other.ptr_;
      other.ptr_ = nullptr;
    }
    return *this;
  }
  intrusive_ptr &operator=(T *ptr) {
    reset(ptr);
    return *this;
  }

  void reset(T *ptr = nullptr) {
    if (ptr != nullptr)
      intrusive_ptr_add_ref(ptr);
    if (ptr_ != nullptr)
      intrusive_ptr_release(ptr_);
    ptr_ = ptr;
  }

  T *get() const { return ptr_; }
  T &operator*() const { return *ptr_; }
  T *operator->() const { return ptr_; }
  explicit operator bool() const { return ptr_ != nullptr; }

private:
  T *ptr_;
};
} // namespace boost

#endif
