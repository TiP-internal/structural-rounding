
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

//Used for ACDS kernel
Set* set_minus(Set*, Set*);
Set* set_union(Set*, Set*);
Set* set_intersection(Set*, Set*);

Set* neighbors(Graph*, int, int);  //pair neighborhood function  
Set* neighbors_closed(Graph*, int);
Set* neighbors_closed(Graph*, int, int);

Set* neighbors_partition1(Graph*, Set*, Set*);
Set* neighbors_partition2(Graph*, Set*, Set*);
Set* neighbors_partition3(Set*, Set*, Set*);


#endif
