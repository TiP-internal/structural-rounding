
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include <vector>

#include "sr_apx/graph/graph.hpp"
#include "sr_apx/setmap/setmap.hpp"
#include "sr_apx/domset/exact/table.hpp"
#include "sr_apx/treewidth/treewidth.hpp"

namespace sr_apx::domset::exact {

enum class Variant {Dom_Set, Indep_Dom_Set, Perf_Dom_Set};

//NOTE driver function, const and opt versions
int calculate(const Graph&, treewidth::Decomposition&, Set&, const Set&, Variant, bool);

void calculate_tables(const Graph&, std::vector<Set>&,
                      std::vector<treewidth::po_bag>&,
                      std::vector<Table*>&,
                      const Set&, Set&, Variant, int);

//---- Calculating Solution
int get_soln_row_index(Table*, const Set&);

//optimization version
int get_solution(std::vector<Table*>&, const Set&);

//constructive version
int get_solution(std::vector<Table*>&, Set&, const Set&);

//constructive and optimization versions
Table* initialize_leaf_table(const Graph&, const Set&, const Set&, treewidth::po_bag, Variant);
Table* join_table(Table*, Table*, const Set&, treewidth::po_bag, bool);
Table* intro_table(const Graph&, Table*, const Set&, const Set&, treewidth::po_bag, Variant, int, bool);
Table* forget_table(Table*, const Set&, treewidth::po_bag, Variant, int, bool);

//-----Helpers
Table* init_parent_table(Table*, bool);

//Dominating Set
int locally_valid_coloring(const Graph&, Table*, const Set&, int, Variant);
int best_intro_child_row(const Set&, std::vector<int>&, const std::vector<int>&, int);
int* best_join_child_rows(Table*, Table*, const Set&, std::vector<int> &);
int get_num_dominators(const Graph&, std::vector<int>&, const std::vector<int>&, int);

}

#endif
