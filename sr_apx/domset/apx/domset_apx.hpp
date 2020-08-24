
#ifndef DOMSET_APX_H
#define DOMSET_APX_H

#include "setmap.hpp"
#include "graph.hpp"

bool is_domset(Graph*, Set*);
int max_deg_vertex(Graph*, Set*);

Set* logn_apx(Graph*);


#endif
