
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include <vector>
#include <iostream>

#include "graph.hpp"
#include "setmap.hpp"
#include "table.hpp"
#include "treewidth.hpp"


enum class Variant {Dom_Set, Indep_Dom_Set, Perf_Dom_Set};

//NOTE driver function, const and opt versions
int calculate(Graph*, TreeDecomp*, Set*, Set*, Variant, bool);  

void calculate_tables(Graph*, std::vector<Set*>&, 
                      std::vector<po_bag>&, 
                      std::vector<Table*>&,
                      Set*, Set*, Variant, int);

//---- Calculating Solution 
int get_soln_row_index(Table*, Set*);

//optimization version
int get_solution(std::vector<Table*>&, Set*);    

//constructive version
int get_solution(std::vector<Table*>&, Set*, Set*);        

//constructive and optimization versions
Table* initialize_leaf_table(Graph*, Set*, Set*, po_bag, Variant);
Table* join_table(Table*, Table*, Set*, po_bag, bool);
Table* intro_table(Graph*, Table*, Set*, Set*, po_bag, Variant, int, bool);
Table* forget_table(Table*, Set*, po_bag, Variant, int, bool);

//-----Helpers
Table* init_parent_table(Table*, bool);

//Dominating Set 
int locally_valid_coloring(Graph*, Table*, Set*, int, Variant);
int best_intro_child_row(Set*, std::vector<int>&, const std::vector<int>&, int);
int* best_join_child_rows(Table*, Table*, Set*, std::vector<int> &);
int get_num_dominators(Graph*, std::vector<int>&, const std::vector<int>&, int);


#endif
