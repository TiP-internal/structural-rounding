
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include <vector> 
#include <deque> 

#include "graph.hpp"
#include "setmap.hpp"
#include "table.hpp"
#include "treedecomp.hpp"


//Normal dominating set is default.
enum class Variant {Dom_Set, Indep_Dom_Set, Perf_Dom_Set};  


Table* calculate_tables(Graph*, std::vector<Set*>&, std::vector<po_bag>&, Set*, Variant);
std::vector<Set*> calc_domset(Graph*, TreeDecomp*, Set*, Variant);

Table* initialize_leaf_table(Graph*, Set*, Set*, Variant);
void update_introduce_table(Graph*, Table*, Set*, Set*, int, Variant);
void update_forget_table(Table*, Set*, int, Variant);
void update_join_table(Table*, Table*, Set*, Variant);
    

//Helpers
void intro_vert_indomset_update(Graph*, Table*, Set*, Row*, int, Variant);
void intro_vert_dominated_update(Graph*, Table*, Set*, Set*, 
                                 Row*, Row*, int, Variant);

int locally_valid_coloring(Graph*, Set*, Row*, std::vector<int> &, Variant);

bool intro_indep_check(Graph*, Table*, std::vector<int>&, int);
bool check_independent(Graph*, Set*);
int get_num_dominators(Graph*, Row*, std::vector<int>&, int);

void minAi_c(Table*, Table*, Set*, Row*, Row*);
int phi(Row*, Set*, std::vector<int>, int);


//For testing
void print_row(Row*);
void print_table(Table*, std::string);
void print_tables(std::vector<Table*>);

#endif
