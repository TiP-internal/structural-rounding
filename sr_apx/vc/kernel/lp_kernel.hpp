
#ifndef LP_KERNEL_H
#define LP_KERNEL_H

#include <tuple>

#include "sr_apx/setmap/setmap.hpp"
#include "sr_apx/graph/graph.hpp"

namespace sr_apx {
namespace vc {
namespace kernel {

std::tuple<Set, Set> lp_kernel(const Graph&);

}}}

#endif
