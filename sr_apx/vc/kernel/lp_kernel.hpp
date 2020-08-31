
#ifndef LP_KERNEL_H
#define LP_KERNEL_H

#include "sr_apx/setmap/setmap.hpp"
#include "sr_apx/graph/graph.hpp"

namespace sr_apx::vc::kernel {

Set** lp_kernel(Graph*);

}

#endif
