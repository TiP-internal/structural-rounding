
#ifndef TREEDECOMP_H
#define TREEDECOMP_H

#include <vector> 
#include <deque> 

#include "graph.hpp"
#include "setmap.hpp"
#include "treewidth.hpp"

class TreeDecomp {
private:
    int tw; 
    int* pre_order;
    
public:
    std::vector<Set*> bags;
    
    std::deque<int> root_indices;  //NOTE: for testing
    std::deque<std::deque<int>> preorder_stack; 
    
    TreeDecomp();
    ~TreeDecomp();
    
    int treewidth();
    int* get_pre_order();
};

#endif
