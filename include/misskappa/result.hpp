#ifndef MISSKAPPA_RESULT_HPP
#define MISSKAPPA_RESULT_HPP

#include <optional>
#include <utility>
#include <variant>

namespace misskappa {

enum class Error {
  invalid_argument,
  dimension_mismatch,
  singular_weight,
  numerical_error,
  not_supported,
  not_converged,
};

// Tag returned by unexpected(); converts implicitly into any Result as its
// error state. Stand-in for std::unexpected<Error> so we stay on C++17.
struct Unexpected {
  Error error;
};

inline Unexpected unexpected(Error e) { return Unexpected{e}; }

// Minimal std::expected<T, Error> stand-in covering exactly what the library
// uses: implicit construction from a value or from unexpected(), a boolean
// has-value test, value access via *, ->, and value(), and error(). Hand-rolled
// (std::expected is C++23) so the package builds on the widest toolchain range.
template <typename T>
class Result {
 public:
  // Value state holding a default-constructed T, matching std::expected (used
  // by the R glue as `Result<T> r;` before assigning the real result).
  Result() = default;
  Result(const T& value) : data_(value) {}
  Result(T&& value) : data_(std::move(value)) {}
  Result(Unexpected u) : data_(u.error) {}

  bool has_value() const { return std::holds_alternative<T>(data_); }
  explicit operator bool() const { return has_value(); }

  T& value() { return std::get<T>(data_); }
  const T& value() const { return std::get<T>(data_); }

  T& operator*() { return std::get<T>(data_); }
  const T& operator*() const { return std::get<T>(data_); }

  T* operator->() { return &std::get<T>(data_); }
  const T* operator->() const { return &std::get<T>(data_); }

  Error error() const { return std::get<Error>(data_); }

 private:
  std::variant<T, Error> data_;
};

// Void specialisation: carries only success/failure plus the error code.
template <>
class Result<void> {
 public:
  Result() = default;
  Result(Unexpected u) : error_(u.error) {}

  bool has_value() const { return !error_.has_value(); }
  explicit operator bool() const { return has_value(); }

  Error error() const { return *error_; }

 private:
  std::optional<Error> error_;
};

}  // namespace misskappa

#endif  // MISSKAPPA_RESULT_HPP
