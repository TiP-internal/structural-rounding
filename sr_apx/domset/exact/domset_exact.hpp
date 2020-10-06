
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include <vector>

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
void add_to_solution(Set*, Row*, std::vector<int>&, Set*);

//optimization version
int get_solution(std::vector<Table*>&, Set*);    

//constructive version
int get_solution(std::vector<Table*>&, Set*);        

//constructive and optimization versions
Table* initialize_leaf_table(Graph*, Set*, Set*, po_bag, Variant);
Table* join_table(Table*, Table*, Set*, po_bag, bool);
Table* intro_table(Graph*, Table*, Set*, Set*, po_bag, Variant, int, bool);
Table* forget_table(Table*, Set*, po_bag, Variant, int, bool);


//-----Helpers
Table* init_parent_table(Table*, bool);
int get_left_join_child_tabind(Set*, std::vector<Table*>&, int);

//Dominating Set 
int locally_valid_coloring(Graph*, Set*, Row*, std::vector<int>&, Variant);
int phi(Row*, Set*, std::vector<int>&, int);
int* minAi_c(Table*, Table*, Set*, Row*);
void intro_vert_indomset_update(Graph*, Table*, Set*, Row*, int, Variant);
int flip_coloring(Table*, std::vector<int>&, std::vector<int>&);
void intro_vert_dominated_update(Graph*, Table*, Set*, Set*, Row*, Row*, int, Variant);

//Independent dominating set + Perfect dominating set
int get_num_dominators(Graph*, Row*, std::vector<int>&, int);
bool intro_indep_check(Graph*, std::vector<int>&, std::vector<int>&, int);
bool check_independent(Graph*, Set*);


#endif
