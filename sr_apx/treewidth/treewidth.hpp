
#ifndef TREEWIDTH_H
#define TREEWIDTH_H

#include <vector>

#include "graph.hpp"
#include "setmap.hpp"


int* treewidth_nodeedit(Graph*, Set*, std::vector<Set*>&, int);
int tree_decomp(Graph*, Set*, Set*);

Set* balanced_separators(Graph*, int);  //greedy alg. from (Althoby et al. 2020)
Set* balanced_separators(Graph*, Set*, int);  //bal. seps for set W.

std::vector<Set*> connected_components(Graph*, Set*);
void dfs(Graph*, Set*, Set*, Set*, int);

int* post_order();

// temporary function to test that vertices from balanced_separators() 
// actually separates graph.
bool test_separators(Graph*, Set*, Set*); 

#endif
