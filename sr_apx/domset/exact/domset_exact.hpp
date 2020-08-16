
#ifndef DOMSET_EXACT_H
#define DOMSET_EXACT_H

#include <vector> 
#include <deque> 

#include "graph.hpp"
#include "setmap.hpp"
#include "table.hpp"
#include "treedecomp.hpp"


Table* calculate_tables(Graph*, std::vector<Set*>&, std::vector<po_bag>&);
std::vector<Set*> calc_domset(Graph*, TreeDecomp*, bool);

//For testing
bool is_domset(Graph*, std::vector<int>);
bool is_ann_domset(Graph*, Set*, Set*);

void print_row(Row*);
void print_table(Table*, std::string);
void print_tables(std::vector<Table*>);
void print_lookups(Table*);

#endif
