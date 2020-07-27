
#ifndef TREEDECOMP_H
#define TREEDECOMP_H

#include <vector> 
#include <deque> 

#include "graph.hpp"
#include "setmap.hpp"
// #include "treewidth.hpp"

class TreeDecomp {
private:
    int tw; 
    std::vector<std::vector<int>> post_order;
    
    void post_order_helper(std::deque<int>&, std::deque<std::deque<int>>&, std::vector<int>&, int);
    bool is_parent(int, std::deque<int> );
    int get_par_index(std::deque<int>, int);
    
public:
    //vector of vectors for each component 
    std::vector<std::vector<Set*>> components_bags;  
    std::vector<std::deque<std::deque<int>>> components_po_stacks; 
    
    TreeDecomp();
    ~TreeDecomp();
    
    //adds roots deque to post_order_stack, and a bags vec for each component
    void add_components(int);  
    
    int treewidth();
    std::vector<std::vector<int>> get_post_order();
    
};

#endif
