#ifndef MISSKAPPA_CORE_TYPES_H
#define MISSKAPPA_CORE_TYPES_H

#include <RcppArmadillo.h>

namespace misskappa_core {

using uvec = arma::Col<arma::uword>;
constexpr int kNaInteger = -2147483648;

} // namespace misskappa_core

#endif // MISSKAPPA_CORE_TYPES_H
