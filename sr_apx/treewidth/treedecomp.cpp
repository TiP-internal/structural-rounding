
#include "treedecomp.hpp"
#include "graph.hpp"

#include <iostream>
#include <fstream>


TreeDecomp::TreeDecomp() {
    tw = -1;
    pre_order = NULL;
    
    std::deque<int> roots;  
    preorder_stack.push_back(roots);
}

TreeDecomp::~TreeDecomp() {
    for(auto ib=bags.begin(); ib!=bags.end(); ib++) {
        delete *ib;
    }
    delete (int*) pre_order;
}

int TreeDecomp::treewidth() {    
    //NOTE: split into get and set methods?
    tw = find_treewidth(bags);          //in treewidth.hpp
    return tw;
}

int* TreeDecomp::get_pre_order() {
    pre_order = new int[bags.size()];
    
    std::deque<int> root_indices = preorder_stack[0]; 
    preorder_stack.pop_front();
    
    int i = 0;
    int rootval = bags.size()-1;  //root bag indices
    
    for(auto iv=root_indices.rbegin(); iv!=root_indices.rend(); iv++) {
        pre_order[i] = rootval;  //visits root first
        
        i++;
        rootval--;
        
        std::deque<int> leaves = preorder_stack.back();
        preorder_stack.pop_back();
        for(auto il=leaves.begin(); il!=leaves.end(); il++) {       
            pre_order[i] = *il;  //then leaves
            i++;
        }
    }
    
    return pre_order;
}

