
#ifndef TREEWIDTH_H
#define TREEWIDTH_H

#include "sr_apx/graph/graph.hpp"
#include "sr_apx/setmap/setmap.hpp"

#include <vector>

namespace sr_apx {
namespace treewidth {

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

public:
    std::vector<Set> components_bags;

    Decomposition();
	Decomposition(const Graph&);

	void build_decomposition(const Graph&);

    int treewidth();
    std::vector<po_bag> get_post_order();
};

Set vertex_delete(const Graph&, int);
Set balanced_separator(const Graph&, const Set&);

}}

#endif
