
#ifndef CDS_PLANAR_KERNEL
#define CDS_PLANAR_KERNEL

#include "setmap.hpp"
#include "graph.hpp"

Graph* planar_kern(Graph*);
Graph* rule1(Graph*);
Graph* rule2(Graph*);
Graph* rule3(Graph*);
Graph* rule4(Graph*);
Graph* rule5(Graph*);

Set* neighbors(Graph*, int, int);  //pair neighborhood function  
Set* neighbors_closed(Graph*, int);
Set* neighbors_closed(Graph*, int, int);

Set* crust(Graph*, Set*, Set*);
Set* mantle(Graph*, Set*, Set*);
Set* core(Set*, Set*, Set*);


#endif
