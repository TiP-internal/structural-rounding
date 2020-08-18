
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include <vector> 
#include <deque> 

#include "graph.hpp"
#include "setmap.hpp"
#include "table.hpp"
#include "treedecomp.hpp"


Table* calculate_tables(Graph*, std::vector<Set*>&, std::vector<po_bag>&, Set*, bool);
std::vector<Set*> calc_domset(Graph*, TreeDecomp*, Set*, bool);

Table* initialize_leaf_table(Graph*, Set*, Set*, int, bool);
void update_introduce_table(Graph*, Table*, Set*, Set*, int, int, bool);
void update_forget_table(Table*, Set*, int, int, bool);
void update_join_table(Table*, Table*, Set*, int, bool);
    

int locally_valid_coloring(Graph*, Set*, Row*, std::vector<int> &, bool);
void minAi_c(Table*, Table*, Set*, Row*, Row*);
int phi(Row*, int);
    

//For testing
bool is_domset(Graph*, std::vector<int>);
bool is_ann_domset(Graph*, Set*, Set*);

void print_row(Row*);
void print_table(Table*, std::string);
void print_tables(std::vector<Table*>);
void print_lookups(Table*);

#endif
