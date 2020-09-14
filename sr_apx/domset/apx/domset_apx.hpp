
#ifndef DOMSET_APX_H
#define DOMSET_APX_H

#include "setmap.hpp"
#include "graph.hpp"

#include <vector>


Set* logn_apx(Graph*);                      // basic greedy alg.
void mod_exp_c_apx(Graph*, Set*, Set*, Set*, int, int);            // [Bourgeois et al. 2009]
Set* Hk_minus_half_apx(Graph*, Set*, Set*, int);       // [Duh & Furer 1997]


//helpers
void semi_local_opt(Graph*, Set*, Set*, Set*, Set*, Set*);
void restricted_phase(Graph*, Set*, Set*, Set*, Set*, Set*, int, int);
int number_one_sets(Graph*, Set*, Set*, Set*);
void maximal_jsets(Graph*, Set*, Set*, Set*, Set*, Set*, int, int);
int vertex_degree(Graph*, Set*, Set*, int);
int vertex_degree(Graph*, Set*, int);
void addto_set(Set*, Set*, Set*);
Set* addto_set(Set*, std::vector<int>, bool);
Set* union_Ssets(Graph*, Set*, Set*, std::vector<int>);

void vertex_combination(Graph*, Set*, Set*, 
                        std::vector<int> &, 
                        std::vector<int> &,
                        std::vector<int> &, 
                        std::vector<int> &, 
                        int, int);

std::vector<int> max_card_sets(Graph*, Set*, Set*, int);
bool p_cardinality(Graph*, Set*, Set*, int);
void single_subset(Graph*, Set*, Set*, Set*);
void included_sets(Graph*, Set*, Set*);

int max_deg_vertex(Graph*, Set*);
int largest_harmonic_num(int);


//for TESTING
bool is_domset(Graph*, Set*);

#endif
