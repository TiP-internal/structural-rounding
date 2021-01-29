
#ifndef VC_LIFT_H
#define VC_LIFT_H

#include "sr_apx/graph/graph.hpp"
#include "sr_apx/setmap/setmap.hpp"

namespace sr_apx {
namespace vc {
namespace lift {

Set naive_lift(const Graph&, const Set&, const Set&);
Set greedy_lift(const Graph&, const Set&, const Set&);
Set apx_lift(const Graph&, const Set&, const Set&);
Set oct_lift(const Graph&, const Set&, const Set&);
Set bip_lift(const Graph&, const Set&, const Set&);
Set recursive_lift(const Graph&, const Set&, const Set&);
Set recursive_oct_lift(const Graph&, const Set&, const Set&);
Set recursive_bip_lift(const Graph&, const Set&, const Set&);

}}}

#endif
