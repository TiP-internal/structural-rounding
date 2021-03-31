
#ifndef DOMSET_LIFT_H
#define DOMSET_LIFT_H

#include "sr_apx/setmap/setmap.hpp"
#include "sr_apx/graph/graph.hpp"

namespace sr_apx {
namespace domset {
namespace lift {

Set naive_lift(const Graph&, const Set&, const Set&);
Set greedy_lift(const Graph&, const Set&, const Set&);

}}}

#endif
