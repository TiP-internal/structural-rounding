
#ifndef TREEDECOMP_H
#define TREEDECOMP_H

#include <vector>
#include <deque>

#include "graph.hpp"
#include "setmap.hpp"
// #include "treewidth.hpp"

struct po_bag {  //for postorder bag.
    int bag_index;
    int num_children;  //leaf if num_children==0
    int parent_bag_index;
    int current_join_child;
};

class TreeDecomp {
private:
    int tw;
    std::vector<std::vector<po_bag>> pre_order;

public:
    //vector of vectors for each component
    std::vector<std::vector<Set*>> components_bags;

    int add_bag(int component, int parent, bool last_child, Set* bag);

    TreeDecomp();
    ~TreeDecomp();

    int treewidth();
    std::vector<std::vector<po_bag>> get_post_order();
};

#endif
