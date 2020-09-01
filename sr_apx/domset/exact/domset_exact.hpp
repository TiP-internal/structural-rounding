
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include <vector> 
#include <deque> 

#include "graph.hpp"
#include "setmap.hpp"
#include "table.hpp"
#include "treedecomp.hpp"


enum class Variant {Dom_Set, Indep_Dom_Set, Perf_Dom_Set};  

//---- Calculating Solution 
int get_soln_row_index(Table*);
void add_to_solution(Set*, Row*, std::vector<int>&);
int get_solution(Table* table);                             //optimization version
Set* get_solution(std::vector<Table*>&, Set*, Row*);        //constructive version

Set* construct_domset(Graph*, TreeDecomp*, Set*, Variant);  //constructive version
int calc_min_domset(Graph*, TreeDecomp*, Set*, Variant);    //optimization version


//--- Constructive Version
Set* treedecomp_reduction(Graph*, std::vector<Set*>&, std::vector<po_bag>);

void calculate_tables(Graph*, std::vector<Set*>&, 
                      std::vector<po_bag>&, std::vector<Table*>&,
                      Set*, Set*, Variant);

Table* initialize_leaf_table(Graph*, Set*, Set*, Variant);
Table* update_introduce_table(Graph*, Table*, Set*, Set*, int, Variant);
Table* update_forget_table(Table*, Set*, int, Variant);
Table* update_join_table(Table*, Table*, Set*, Variant);


//----- Updates the child table to be par table.
void merge_introduce_table(Graph*, Table*, Set*, Set*, int, Variant);
void merge_forget_table(Table*, Set*, int, Variant);
void merge_join_table(Table*, Table*, Set*, Variant);


//Helpers
int locally_valid_coloring(Graph*, Set*, Row*, std::vector<int>&, Variant);
int get_num_dominators(Graph*, Row*, std::vector<int>&, int);
bool intro_indep_check(Graph*, std::vector<int>&, std::vector<int>&, int);
bool check_independent(Graph*, Set*);
int phi(Row*, Set*, std::vector<int>&, int);
void remove_node_from_postack(std::vector<po_bag> &, po_bag &); 
void remove_edge_from_postack(std::vector<po_bag> &, po_bag &);
bool is_special_subset(Set*, Set*, Set*);
int* minAi_c(Table*, Table*, Set*, Row*, Row*);
void intro_vert_indomset_update(Graph*, Table*, Set*, Row*, int, Variant);
void intro_vert_dominated_update(Graph*, Table*, Set*, Set*, Row*, Row*, int, Variant);

//-----For testing
void print_row(Row*);
void print_table(Table*, std::string);

void print_tables(std::vector<Table*>&);
void print_postorder(std::vector<po_bag>);
void print_pobag(po_bag);

#endif
