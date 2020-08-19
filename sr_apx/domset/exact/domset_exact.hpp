
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include <vector> 
#include <deque> 

#include "graph.hpp"
#include "setmap.hpp"
#include "table.hpp"
#include "treedecomp.hpp"


Table* calculate_tables(Graph*, std::vector<Set*>&, std::vector<po_bag>&, Set*);
int calc_domset(Graph*, TreeDecomp*, Set*);

Table* initialize_leaf_table(Graph*, Set*, Set*);
void update_introduce_table(Graph*, Table*, Set*, Set*, int);
void update_forget_table(Table*, Set*, int);
void update_join_table(Table*, Table*, Set*);
    

int locally_valid_coloring(Graph*, Set*, Row*, std::vector<int> &);
void minAi_c(Table*, Table*, Set*, Row*, Row*);
int phi(Row*, Set*, std::vector<int>, int);
    

//For testing
bool is_domset(Graph*, std::vector<int>);
bool is_ann_domset(Graph*, Set*, Set*);

void print_row(Row*);
void print_table(Table*, std::string);
void print_tables(std::vector<Table*>);
//void print_lookups(Table*);

#endif
