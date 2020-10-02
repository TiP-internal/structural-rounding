
#include "treedecomp.hpp"
#include "graph.hpp"

#include <iostream>
#include <fstream>


TreeDecomp::TreeDecomp() {
    tw = -1;
}

TreeDecomp::~TreeDecomp() {
    for (Set* s : components_bags) {
        delete s;
    }
}

int TreeDecomp::treewidth() {
    return tw-1;
}

std::vector<po_bag> TreeDecomp::get_post_order() {
    std::vector<po_bag> po;
    for (auto it = pre_order.rbegin(); it != pre_order.rend(); ++it) {
        po.push_back(*it);
    }

    return po;
}

Set* copyset(Set* s) {
    Set* copy = new Set();
    for (int x : *s) {
        copy->insert(x);
    }
    return copy;
}

// returns placement of new bag
// parent_index should be 0 for root bag
// only modifies pre_order and components_bags
int TreeDecomp::add_bag(int parent, bool last_child, Set* bag) {
    // ensures that root bag exists
    if (components_bags.empty()) {
        components_bags.push_back(new Set());
        po_bag p;
        p.bag_index = 0;
        p.num_children = 0;
        p.parent_bag_index = -1;
        p.current_join_child = 0;
        pre_order.push_back(p);
    }

    int true_parent = pre_order[parent].current_join_child;

    if (pre_order[true_parent].num_children > 0 && !last_child) {
        // create a new join bag to hold the child
        int index = components_bags.size();
        components_bags.push_back(copyset(components_bags[true_parent]));
        po_bag p;
        p.bag_index = index;
        p.num_children = 0;
        p.parent_bag_index = true_parent;
        p.current_join_child = index;
        pre_order.push_back(p);
        pre_order[true_parent].num_children++;
        pre_order[parent].current_join_child = index;
        true_parent = index;
    }

    if (pre_order[true_parent].num_children > 0 || !last_child) {
        int index = components_bags.size();
        components_bags.push_back(copyset(components_bags[true_parent]));
        po_bag p;
        p.bag_index = index;
        p.num_children = 0;
        p.parent_bag_index = true_parent;
        p.current_join_child = index;
        pre_order.push_back(p);
        pre_order[true_parent].num_children++;
        true_parent = index;
    }

    // create the child
    int use_parent = true_parent;
    Set* use_bag = copyset(components_bags[true_parent]);

    // add introduce bags
    for (int x : *(components_bags[true_parent])) {
        // only look at vertices to be added
        if (bag->contains(x)) {
            continue;
        }

        use_bag->erase(x);
        int index = components_bags.size();
        components_bags.push_back(copyset(use_bag));
        po_bag p;
        p.bag_index = index;
        p.num_children = 0;
        p.parent_bag_index = use_parent;
        p.current_join_child = index;
        pre_order.push_back(p);
        pre_order[use_parent].num_children++;
        use_parent = index;
    }

    // add forget bags
    for (int x : *bag) {
        // only look at vertices to be removed
        if (components_bags[true_parent]->contains(x)) {
            continue;
        }

        use_bag->insert(x);
        int index = components_bags.size();
        components_bags.push_back(copyset(use_bag));
        po_bag p;
        p.bag_index = index;
        p.num_children = 0;
        p.parent_bag_index = use_parent;
        p.current_join_child = index;
        pre_order.push_back(p);
        pre_order[use_parent].num_children++;
        use_parent = index;
    }

    tw = tw < bag->size() ? bag->size() : tw;

    // return index of final bag
    return components_bags.size() - 1;
}
