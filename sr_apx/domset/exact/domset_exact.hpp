
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
template<class T>
int get_soln_row_index(Table<T>*);
void add_to_solution(Set*, RowConstruct*, std::vector<int>&);
int get_solution(Table<Row*>* table);                                       //decision version
Set* get_solution(std::vector<Table<RowConstruct*>*>&, Set*, RowConstruct*);//constructive version

Set* construct_domset(Graph*, TreeDecomp*, Set*, Variant);  //constructive version
int calc_min_domset(Graph*, TreeDecomp*, Set*, Variant);    //decision version


//--- Constructive Version
void calculate_tables(Graph*, std::vector<Set*>&, 
                      std::vector<po_bag>&, std::vector<Table<RowConstruct*>*>&,
                      Set*, Variant);

Table<RowConstruct*>* initialize_leaf_table_const(Graph*, Set*, Set*, Variant);

Table<RowConstruct*>* update_introduce_table(Graph*, Table<RowConstruct*>*,
                                             Set*, Set*, int, Variant);

Table<RowConstruct*>* update_forget_table(Table<RowConstruct*>*, Set*, int, Variant);

Table<RowConstruct*>* update_join_table(Table<RowConstruct*>*, 
                                        Table<RowConstruct*>*, Set*, Variant);


//----- Decision version
Table<Row*>* calculate_tables(Graph*, std::vector<Set*>&, std::vector<po_bag>&, Set*, Variant); 
Table<Row*>* initialize_leaf_table(Graph*, Set*, Set*, Variant);
void update_introduce_table(Graph*, Table<Row*>*, Set*, Set*, int, Variant);
void update_forget_table(Table<Row*>*, Set*, int, Variant);
void update_join_table(Table<Row*>*, Table<Row*>*, Set*, Variant);
    


//Helpers
int locally_valid_coloring(Graph*, Set*, Row*, std::vector<int>&, Variant);
int get_num_dominators(Graph*, Row*, std::vector<int>&, int);
bool intro_indep_check(Graph*, std::vector<int>&, std::vector<int>&, int);
bool check_independent(Graph*, Set*);
int phi(Row*, Set*, std::vector<int>&, int);

template<class T>
int* minAi_c(Table<T>*, Table<T>*, Set*, Row*, Row*);

template<class T>
void intro_vert_indomset_update(Graph*, Table<T>*, Set*, T, int, Variant);

template<class T>
void intro_vert_dominated_update(Graph*, Table<T>*, Set*, Set*, 
                                 T, T, int, Variant);


//-----For testing
void print_row(Row*);
void print_row(RowConstruct*);

template<class T>
void print_table(Table<T>*, std::string);

template<class T>
void print_tables(std::vector<Table<T>*>&);

#endif
