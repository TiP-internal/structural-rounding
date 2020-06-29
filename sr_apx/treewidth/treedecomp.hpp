
#ifndef TREEDECOMP_H
#define TREEDECOMP_H

#include <vector> 

#include "graph.hpp"
#include "setmap.hpp"
#include "treewidth.hpp"

class TreeDecomp {
private:
    int tw; 
    int* post_order;
    
public:
    std::vector<Set*> bags;
    
    TreeDecomp();
    ~TreeDecomp();
    
    int treewidth();
    int* get_post_order();
};

#endif
