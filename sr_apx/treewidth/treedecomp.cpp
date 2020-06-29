
#include "treedecomp.hpp"
#include "graph.hpp"

#include <iostream>
#include <fstream>


TreeDecomp::TreeDecomp(Graph* g, Set* V) {
    tw = -1;
    graph = g;
    post_order = NULL;
    
    Z = V;  //??
    W = new Set();
}

TreeDecomp::~TreeDecomp() {
    for(auto ib=bags.begin(); ib!=bags.end(); ib++) {
        delete *ib;
    }
    delete Z, W;
    delete (int*) post_order;
}

void TreeDecomp::tree_decomposition() {
    tree_decomp(graph, Z, W, bags);  //in treewidth.hpp
}

int TreeDecomp::treewidth() {    
    tw = find_treewidth(bags);
    return tw;
}

int* TreeDecomp::get_post_order() {
    printf("returning post_order()\n");
    return 0;
}

