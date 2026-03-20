#ifndef _PUBLIC_PERMUTATION_HPP
#define _PUBLIC_PERMUTATION_HPP

#include "orbit/internal/move/permutation_impl.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_permutation.hpp"

namespace orbit {

using permutation = permutation_impl<>;

} // namespace orbit

#endif