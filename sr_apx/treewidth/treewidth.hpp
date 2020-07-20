
#ifndef TREEWIDTH_H
#define TREEWIDTH_H

#include <vector>
#include <deque> 

#include "graph.hpp"
#include "setmap.hpp"

Set* treewidth_nodeedit(Graph*, int);
Set* tree_decomp(Graph*, Set*, Set*, std::deque<std::deque<int>>&, std::vector<Set*>&);
int find_treewidth(std::vector<Set*> &);

Set* balanced_separators(Graph*, int);        //greedy alg. from (Althoby et al. 2020)
Set* balanced_separators(Graph*, Set*, int);  //bal. seps for set W.

std::vector<Set*> connected_components(Graph*);
void dfs(Graph*, Set*, Set*, int);
int min_deg_vert(Graph*);
bool sets_equal(Set*, Set*);

bool test_separators(Graph*, Set*, Set*, int, int, int);
#endif
