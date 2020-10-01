
#include "treedecomp.hpp"
#include "graph.hpp"

#include <iostream>
#include <fstream>


TreeDecomp::TreeDecomp() {
    tw = -1;
    nice = false;
}

TreeDecomp::~TreeDecomp() {
    for(auto ib=components_bags.begin(); ib!=components_bags.end(); ib++) {
        std::vector<Set*> comp = *ib;
        for(auto ic=comp.begin(); ic!=comp.end(); ic++) {
            delete *ic;
        }
    }
}

void TreeDecomp::add_components(int n_components) {
    /*
     * Adds a component to the decomposition.
     * Need separate stacks for each component.
     */

    for(int i=0; i<n_components; i++) {
        std::deque<int> roots;
        std::deque<std::deque<int>> po_comp_stack;  //first deque are the roots of decomp of component
        po_comp_stack.push_back(roots);
        components_po_stacks.push_back(po_comp_stack);

        std::vector<Set*> bags_comp;
        components_bags.push_back(bags_comp);
    }
}

int TreeDecomp::treewidth() {
    int tw = 0;
    for(auto ib=components_bags.begin(); ib!=components_bags.end(); ib++) {  //for each comp.
        std::vector<Set*> bag_comp = *ib;
        for(auto ic=bag_comp.begin(); ic!=bag_comp.end(); ic++) {
            Set* bag = *ic;
            if (bag->size() > tw) tw = bag->size();
        }
    }
    return tw-1;
}

bool TreeDecomp::is_parent(int leaf, std::deque<int> pars) {
    for(auto ip=pars.begin(); ip!=pars.end(); ip++) {
        int par = *ip;
        if(leaf == par) return true;
    }
    return false;
}

int TreeDecomp::get_par_index(std::deque<int> pars, int par_node) {
    int index = -1;
    for(int j=0; j<pars.size(); j++) {
        if(pars[j] == par_node) index = j;
    }
    return index;
}

std::vector<std::vector<po_bag>> TreeDecomp::get_post_order() {
    if (nice) {
        std::vector<std::vector<po_bag>> po;
        for (std::vector<po_bag> comp : post_order) {
            std::vector<po_bag> rev;
            for (auto it = comp.rbegin(); it != comp.rend(); ++it) {
                rev.push_back(*it);
                rev[rev.size() - 1].parent_bag_index = comp.size() - 1 - (*it).parent_bag_index;
            }
            po.push_back(rev);
        }

        post_order = po;
        return post_order;
    }

    /*
     * Finds the post_order for each component.
     */
    int i=0;
    for(auto ib=components_bags.begin(); ib!=components_bags.end(); ib++) {  //for each comp.
        std::vector<Set*> bag_comp = *ib;
        std::deque<std::deque<int>> po_stack = components_po_stacks[i];

        std::vector<po_bag> po_comp;
        std::deque<int> pars = po_stack[0];
        po_stack.pop_front();

        post_order_helper(pars, po_stack, po_comp, pars[0], -1);
        post_order.push_back(po_comp);
        i++;
    }

    return post_order;
}

void TreeDecomp::post_order_helper(std::deque<int> &pars,
                                   std::deque<std::deque<int>> &leaf_stack,
                                    std::vector<po_bag>& po,
                                   int parent, int grandparent) {
    /*
     * Finds post order for the given component.
     *
     * Input format:
     *      - pars: deque of the parent nodes
     *          ex: [par0, par1, ...]
     *      - leaf_stack: deque of deques. each deque contains the leaf
     *                    nodes of the parent node at that index.
     *          ex: [ [par0_leaf0, par0_leaf1, ...],
     *                [par1_leaf0, par1_leaf1, ...], ...]
     *
     */
    int parent_ind = parent;
    int par_i = get_par_index(pars, parent_ind);
    //TODO could add this to recursive call

    //printf("parent_ind=%d, par_i=%d\n", parent_ind, par_i);

    std::deque<int> leaves = leaf_stack[par_i];

    //-------left+right  first
    for(auto il=leaves.begin(); il!=leaves.end(); il++) {
        int leaf_ind = *il;
        bool ispar = is_parent(leaf_ind, pars); //is the leaf a parent also?
        if(ispar) {
            post_order_helper(pars, leaf_stack, po, leaf_ind, parent_ind);
        } else {
            po_bag lf;
            lf.num_children=0;
            lf.parent_bag_index = parent_ind;
            lf.bag_index=leaf_ind;
            po.push_back(lf);
        }
    }

    //------roots last
    po_bag prnt;
    prnt.num_children=leaves.size();
    prnt.bag_index=parent_ind;
    prnt.parent_bag_index=grandparent;
    po.push_back(prnt);
}

void TreeDecomp::set_nice() {
    nice = true;
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
// only modifies post_order and components_bags
int TreeDecomp::add_bag(int component, int parent, bool last_child, Set* bag) {
    // ensures that component exists
    if (component >= components_bags.size()) {
        components_bags.resize(component + 1);
        post_order.resize(component + 1);
    }

    // ensures that empty root bag is placed
    if (components_bags[component].empty()) {
        components_bags[component].push_back(new Set());
        po_bag p;
        p.bag_index = 0;
        p.num_children = 0;
        p.parent_bag_index = -1;
        p.current_join_child = 0;
        post_order[component].push_back(p);
    }

    int true_parent = post_order[component][parent].current_join_child;

    if (post_order[component][true_parent].num_children > 0 && !last_child) {
        // create a new join bag to hold the child
        int index = components_bags[component].size();
        components_bags[component].push_back(copyset(components_bags[component][true_parent]));
        po_bag p;
        p.bag_index = index;
        p.num_children = 0;
        p.parent_bag_index = true_parent;
        p.current_join_child = index;
        post_order[component].push_back(p);
        post_order[component][true_parent].num_children++;
        post_order[component][parent].current_join_child = index;
        true_parent = index;
    }

    if (post_order[component][true_parent].num_children > 0 || !last_child) {
        int index = components_bags[component].size();
        components_bags[component].push_back(copyset(components_bags[component][true_parent]));
        po_bag p;
        p.bag_index = index;
        p.num_children = 0;
        p.parent_bag_index = true_parent;
        p.current_join_child = index;
        post_order[component].push_back(p);
        post_order[component][true_parent].num_children++;
        true_parent = index;
    }

    // create the child
    int use_parent = true_parent;
    Set* use_bag = copyset(components_bags[component][true_parent]);

    // add introduce bags
    for (int x : *(components_bags[component][true_parent])) {
        // only look at vertices to be added
        if (bag->contains(x)) {
            continue;
        }

        use_bag->erase(x);
        int index = components_bags[component].size();
        components_bags[component].push_back(copyset(use_bag));
        po_bag p;
        p.bag_index = index;
        p.num_children = 0;
        p.parent_bag_index = use_parent;
        p.current_join_child = index;
        post_order[component].push_back(p);
        post_order[component][use_parent].num_children++;
        use_parent = index;
    }

    // add forget bags
    for (int x : *bag) {
        // only look at vertices to be removed
        if (components_bags[component][true_parent]->contains(x)) {
            continue;
        }

        use_bag->insert(x);
        int index = components_bags[component].size();
        components_bags[component].push_back(copyset(use_bag));
        po_bag p;
        p.bag_index = index;
        p.num_children = 0;
        p.parent_bag_index = use_parent;
        p.current_join_child = index;
        post_order[component].push_back(p);
        post_order[component][use_parent].num_children++;
        use_parent = index;
    }

    // return index of final bag
    return components_bags[component].size() - 1;
}

void TreeDecomp::make_nice_helper(int component, int node, int nice_parent, bool last_child, TreeDecomp* nice_decomp) {
    int x = nice_decomp->add_bag(component, nice_parent, last_child, components_bags[component][node]);

    int index = get_par_index(components_po_stacks[component][0], node);
    if (index == -1) {
        nice_decomp->add_bag(component, x, true, new Set());
        return;
    }

    int size = components_po_stacks[component][index + 1].size();
    for (int i = 0; i < size - 1; ++i) {
        make_nice_helper(component, components_po_stacks[component][index + 1][i], x, false, nice_decomp);
    }

    if (size > 0) {
        int last = components_po_stacks[component][index + 1][size - 1];
        make_nice_helper(component, last, x, true, nice_decomp);
    }
}

TreeDecomp* TreeDecomp::make_nice() {
    TreeDecomp* nice_decomp = new TreeDecomp();
    nice_decomp->nice = true;

    for (int component = 0; component < components_bags.size(); ++component) {
        make_nice_helper(component, components_po_stacks[component][0][0], 0, true, nice_decomp);
    }

    return nice_decomp;
}
