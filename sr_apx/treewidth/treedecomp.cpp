
#include "treedecomp.hpp"
#include "graph.hpp"
#include "util.hpp"

#include <iostream>
#include <fstream>


TreeDecomp::TreeDecomp(Graph* g, int tw) {
    target_tw = tw;
    
    edit_set = NULL;
    graph = g;
    post_order = NULL;
}

TreeDecomp::~TreeDecomp() {
    for(auto ib=bags.begin(); ib!=bags.end(); ib++) {
        delete *ib;
    }
    delete (int*) post_order;
}

Set* TreeDecomp::treewidth_edit() {
    post_order = treewidth_nodeedit(graph, edit_set, bags, target_tw);  //in treewidth.hpp
    
    return edit_set;
}

