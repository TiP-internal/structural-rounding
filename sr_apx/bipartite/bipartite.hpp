
#ifndef BIPARTITE_H
#define BIPARTITE_H

#include "sr_apx/graph/graph.hpp"
#include "sr_apx/setmap/setmap.hpp"

Set* vertex_delete(Graph*);
Set* prescribed_octset(Graph*, const char*);
Set** verify_bipartite(Graph*, Set*);

#endif
