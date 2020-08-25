
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include <vector> 
#include <deque> 

#include "graph.hpp"
#include "setmap.hpp"
#include "table.hpp"
#include "treedecomp.hpp"


//Normal dominating set is default.
enum class Variant {Indep_Dom_Set, Perf_Dom_Set};  


Table* calculate_tables(Graph*, std::vector<Set*>&, std::vector<po_bag>&, Set*, Variant);
std::vector<Set*> calc_domset(Graph*, TreeDecomp*, Set*, Variant);

Table* initialize_leaf_table(Graph*, Set*, Set*, Variant);
void update_introduce_table(Graph*, Table*, Set*, Set*, int, Variant);
void update_forget_table(Table*, Set*, int, Variant);
void update_join_table(Table*, Table*, Set*, Variant);
    

int locally_valid_coloring(Graph*, Set*, Row*, std::vector<int> &, Variant);
bool check_independent(Graph*, Set*);

void minAi_c(Table*, Table*, Set*, Row*, Row*);
int phi(Row*, Set*, std::vector<int>, int);
    

//---For testing
bool is_domset(Graph*, std::vector<int>);
bool is_ann_domset(Graph*, Set*, Set*);
bool is_indepen_domset(Graph*, Set*);
bool is_indepen_ann_domset(Graph*, Set*, Set*); 
bool is_perf_domset(Graph*, Set*);
bool is_per_ann_domset(Graph*, Set*, Set*);

void print_row(Row*);
void print_table(Table*, std::string);
void print_tables(std::vector<Table*>);

#endif
