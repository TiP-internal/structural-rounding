
#ifndef MATCHING_H
#define MATCHING_H

#include "sr_apx/setmap/setmap.hpp"
#include "sr_apx/graph/graph.hpp"

namespace sr_apx::misc {

Map<int>* bipartite_matching(Graph*, Set*, Set*);

}

#endif
