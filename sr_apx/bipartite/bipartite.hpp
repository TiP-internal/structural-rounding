
#ifndef BIPARTITE_H
#define BIPARTITE_H

#include "sr_apx/graph/graph.hpp"
#include "sr_apx/setmap/setmap.hpp"

namespace sr_apx::bipartite {

Set* vertex_delete(Graph*);
Set* prescribed_octset(Graph*, const char*);
Set** verify_bipartite(Graph*, Set*);

}

#endif
