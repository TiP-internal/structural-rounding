
#ifndef TREEWIDTH_H
#define TREEWIDTH_H

#include "sr_apx/graph/graph.hpp"
#include "sr_apx/setmap/setmap.hpp"

#include <vector>

namespace sr_apx {
namespace treewidth {

struct po_bag {
    int left_child = -1;
    int right_child = -1;
    int current_join_child = 0;

    std::vector<int> bag;

    po_bag() = default;
    po_bag(int index, const std::vector<int>& b): current_join_child(index), bag(b) {}
    po_bag(int index, std::vector<int>&& b): current_join_child(index), bag(std::move(b)) {}
};

class Decomposition {
private:
    bool build;
    int width;

    int add_bag(int parent, bool last_child, const Set& bag);

    void tree_decomp(Graph&, Set&, Set&, int, bool);
    void tree_decomp_ordering(Graph&, int, std::vector<int>);


public:
    std::vector<po_bag> pre_order;

    explicit Decomposition(bool b);
    Decomposition(const Graph&);

    void build_decomposition(const Graph&);
    void build_decomposition(const Graph&, int);
    void build_decomposition(const Graph&, std::vector<int>);

    void build_gd_decomposition(const Graph&);
    void build_gd_decomposition(const Graph&, int);
    void build_gfi_decomposition(const Graph&);
    void build_gfi_decomposition(const Graph&, int);

    void sort_bags();

    int treewidth();
};

Set vertex_delete(const Graph&, int);
Set vertex_gd_delete(const Graph&, int);
Set vertex_gfi_delete(const Graph&, int);
Set balanced_separator(const Graph&, const Set&);
Set close_balanced_separator(const Graph&, const Set&);

int fill_edges(const Graph&, int);
int min_vertex(const Graph&, std::vector<Set>&, int);
std::vector<int> greedy_fill_in(const Graph&, int);
std::vector<int> greedy_degree(const Graph&, int);
Graph fill(const Graph&, int, std::vector<int>);
Graph minimal_triangulation(const Graph&, int);

}}

#endif
