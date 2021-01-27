
#ifndef TREEWIDTH_H
#define TREEWIDTH_H

#include "sr_apx/graph/graph.hpp"
#include "sr_apx/setmap/setmap.hpp"

#include <vector>

namespace sr_apx::treewidth {

struct po_bag {  //for postorder bag.
    int bag_index;
    int num_children = 0;  //leaf if num_children==0
    int parent_bag_index;
    int current_join_child;

    po_bag() = default;
    po_bag(int index, int parent): bag_index(index), parent_bag_index(parent), current_join_child(index) {}
};

class Decomposition {
private:
    int tw;
    std::vector<po_bag> pre_order;

	int add_bag(int parent, bool last_child, Set& bag);
	void tree_decomp(Graph&, Set&, Set&, int, bool);
    void tree_decomp_ordering(Graph&, int, std::vector<int>);


public:
    std::vector<Set> components_bags;

    Decomposition();
	Decomposition(const Graph&);

	void build_decomposition(const Graph&);
    void build_decomposition(const Graph&, std::vector<int>);

    int treewidth();
    std::vector<po_bag> get_post_order();
};

Set vertex_delete(const Graph&, int);
Set balanced_separator(const Graph&, const Set&);

std::vector<int> greedy_fill_in(const Graph&, int);
std::vector<int> greedy_degree(const Graph&, int);
Graph fill(const Graph&, int, std::vector<int>);
Graph minimal_triangulation(const Graph&, int);

}

#endif
