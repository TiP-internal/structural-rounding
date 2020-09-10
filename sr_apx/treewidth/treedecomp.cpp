
#include "treedecomp.hpp"
#include "graph.hpp"

#include <iostream>
#include <fstream>


TreeDecomp::TreeDecomp() {
    tw = -1;
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

//NOTE tracks the number of recursive calls made. For each call, we take a new parent from pars.
//The index of this corresponds to the index of leaves in the leaf_stack. 
int NUM_CALLS=0;  
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
    
    //NOTE this is doing the same thing as NUM_CALLS, but takes much more time
    //int par_i = get_par_index(pars, parent_ind); 
    
    if(NUM_CALLS>=leaf_stack.size()) {
        printf("WARNING: NUM_CALLS in TreeDecomp::post_order_helper incorrect.\n");
    } 
    
    std::deque<int> leaves = leaf_stack[NUM_CALLS];
    NUM_CALLS++;
    
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















