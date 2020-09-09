
#ifndef DOMSET_APX_H
#define DOMSET_APX_H

#include "setmap.hpp"
#include "graph.hpp"


Set* logn_apx(Graph*);              // basic greedy alg.
Set* mod_exp_c_apx(Graph*, int);    // [Bourgeois et al. 2009]
Set* Hk_minus_half_apx(Graph*, Set*);     // [Duh & Furer 1997]


//helpers
int max_deg_vertex(Graph*, Set*);
int largest_harmonic_num(int);


//for TESTING
bool is_domset(Graph*, Set*);

#endif
