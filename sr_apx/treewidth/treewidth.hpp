
#ifndef TREEWIDTH_H
#define TREEWIDTH_H

#include "graph.hpp"
#include "setmap.hpp"

#include <vector>

struct po_bag {  //for postorder bag.
    int bag_index;
    int num_children = 0;  //leaf if num_children==0
    int parent_bag_index;
    int current_join_child;

    po_bag() = default;
    po_bag(int index, int parent): bag_index(index), parent_bag_index(parent), current_join_child(index) {}
};

class TreeDecomp {
private:
    int tw;
    std::vector<po_bag> pre_order;

	int add_bag(int parent, bool last_child, Set* bag);
	void tree_decomp(Graph*, Set*, Set*, int, bool);

public:
    //vector of vectors for each component
    std::vector<Set*> components_bags;

    TreeDecomp();
	TreeDecomp(Graph*);
    ~TreeDecomp();

	void build_decomposition(Graph*);

    int treewidth();
    std::vector<po_bag> get_post_order();
};

Set* treewidth_nodeedit(Graph*, Set*, int, bool);

Set* balanced_separators(Graph*, int);        //greedy alg. from (Althoby et al. 2020)
Set* balanced_separators(Graph*, Set*, int);  //bal. seps for set W.

std::vector<Set*> connected_components(Graph*);

int min_deg_vert(Graph*);

//For testing
bool sets_equal(Set*, Set*);
bool test_separators(Graph*, Set*, Set*, int, int, int);
bool is_clique(Graph*);
#endif
