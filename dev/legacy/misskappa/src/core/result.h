#ifndef MISSKAPPA_CORE_RESULT_H
#define MISSKAPPA_CORE_RESULT_H

#include <optional>
#include <string>
#include <utility>

namespace misskappa_core {

enum class Status { kOk, kError };

template <typename T>
struct Result {
  Status status = Status::kError;
  std::optional<T> value;
  std::string error_message;

  bool IsOk() const { return status == Status::kOk; }

  static Result Ok(T v) { return {Status::kOk, std::move(v), ""}; }
  static Result Error(std::string msg) { return {Status::kError, std::nullopt, std::move(msg)}; }
};

} // namespace misskappa_core

#endif // MISSKAPPA_CORE_RESULT_H
