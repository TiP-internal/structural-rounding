
#ifndef TREEDECOMP_H
#define TREEDECOMP_H

#include <vector> 

#include "graph.hpp"
#include "setmap.hpp"
#include "treewidth.hpp"

class TreeDecomp {
public:
    int target_tw;
    
    Graph* graph;
    Set* edit_set;
    
    std::vector<Set*> bags;  
    int* post_order;
    
    TreeDecomp(Graph*, int);
    ~TreeDecomp();
    
    Set* treewidth_edit();
};

#endif
