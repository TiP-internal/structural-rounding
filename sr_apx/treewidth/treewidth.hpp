
#ifndef TREEWIDTH_H
#define TREEWIDTH_H

#include "graph.hpp"
#include "setmap.hpp"


Set* treewidth_nodeedit(Graph*, int);
int tree_decomp(Graph*, Set*, Set*);

Set* balanced_separators(Graph*, int);  //greedy alg. from (Althoby et al. 2020)
Set* balanced_separators(Graph*, Set*, int);  //bal. seps for set W.

// temporary function to test that vertices from balanced_separators() 
// actually separates graph.
bool test_separators(Graph*, Set*, Set*); 

#endif
