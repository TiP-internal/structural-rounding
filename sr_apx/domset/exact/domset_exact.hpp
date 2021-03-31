
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include "sr_apx/graph/graph.hpp"
#include "sr_apx/setmap/setmap.hpp"
#include "sr_apx/treewidth/treewidth.hpp"

namespace sr_apx {
namespace domset {
namespace exact {

Set tw_exact(const Graph&, treewidth::Decomposition&, const Set&);

}}}

#endif
