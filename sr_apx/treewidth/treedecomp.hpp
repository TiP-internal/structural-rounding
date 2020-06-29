
#ifndef TREEDECOMP_H
#define TREEDECOMP_H

#include <vector> 

#include "graph.hpp"
#include "setmap.hpp"
#include "treewidth.hpp"

class TreeDecomp {
public:
    int tw; 
    
    Graph* graph;
    
    Set* Z;
    Set* W;
    
    std::vector<Set*> bags;  
    int* post_order;
    
    TreeDecomp(Graph*, Set*);
    ~TreeDecomp();
    
    void tree_decomposition();
    int treewidth();
    int* get_post_order();
};

#endif
