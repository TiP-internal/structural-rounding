
#ifndef BIPARTITE_H
#define BIPARTITE_H

#include <tuple>

#include "sr_apx/graph/graph.hpp"
#include "sr_apx/setmap/setmap.hpp"

namespace sr_apx::bipartite {

Set vertex_delete(const Graph&);
Set prescribed_octset(const Graph&, const char*);
std::tuple<Set, Set, Set> verify_bipartite(const Graph&, const Set&);

}

#endif
