
#include "treedecomp.hpp"
#include "graph.hpp"

#include <iostream>
#include <fstream>


TreeDecomp::TreeDecomp() {
    tw = -1;
    post_order = NULL;
}

TreeDecomp::~TreeDecomp() {
    for(auto ib=bags.begin(); ib!=bags.end(); ib++) {
        delete *ib;
    }
    delete (int*) post_order;
}

int TreeDecomp::treewidth() {    
    //NOTE: split into get and set methods?
    tw = find_treewidth(bags);          //in treewidth.hpp
    return tw;
}

int* TreeDecomp::get_post_order() {
    post_order = new int[bags.size()];
    
    int index = 0;
    for(int i=bags.size(); i>=0; i--) {
        post_order[index] = i;
        index++;
    }
    return post_order;
}

