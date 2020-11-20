
#ifndef DOMSET_APX_H
#define DOMSET_APX_H

#include "sr_apx/setmap/setmap.hpp"
#include "sr_apx/graph/graph.hpp"

namespace sr_apx {
namespace domset {
namespace apx {

Set greedy_apx(const Graph&);

}}}

// struct Domset{
//     Set* base = new Set();  //universe
//
//     //map of subset of the univers. key=vert ,
//     //Set* are neighbors of vert + vert
//     Map<Set*>* collection = new Map<Set*>();
//
//     Set* domset = new Set();
// };
//
//
// struct OPT_2_1_Sets {
//     Map<Set*>* two_sets_added = new Map<Set*>();
//     Map<Set*>* one_sets_added = new Map<Set*>();
// };
//
//
// //----apx algs
// // basic greedy alg.
// Set* logn_apx(Graph*);
// int logn_apx(Domset&);
//
// // [Bourgeois et al. 2009]
// Set* mod_exponential_domset(Domset&, int);
// int mod_exp_c_apx(Domset, int);
//
// // [Duh & Furer 1997]
// int Hk_minus_half_apx(Domset&, int);
//
//
// //----misc
// Domset reduction(Graph*);
//
//
// //----helpers
// bool is_subset(Set*, Set*);
// void remove_from_collection(Domset&, int);
// int setsize(Map<Set*>*, Set*);
// int how_many(Domset&, std::vector<int>&);
//
// //--hk helpers
// int maximal_jsets(Domset&, int, int);
// int maximal_set(Domset&, int);
//
// int restricted_phase(Domset&, int, int);
// int semi_local_opt(Domset&, int);
//
// void greedy_3set_select(Domset&, Map<Set*>*);
// OPT_2_1_Sets opt_2_1_set_cover(Domset&, Map<Set*>*);
// int st_improvements(Domset&, Map<Set*>*, const int);
// std::vector<std::vector<int>> threeset_combs(Map<Set*>*);
//
// void max_matching(Graph*, Domset&, OPT_2_1_Sets&, Map<Set*>*);
// bool does_del_incr_onesets(Domset&, int);
//
//
// //---mod exp helpers
// Domset remove_sets_from_S(Domset, std::vector<int>);
// Domset remove_from_SandC_addsoln(Domset, std::vector<int>);
//
// void vertex_combination(Domset&,
//                         std::vector<int> &,
//                         std::vector<int> &,
//                         std::vector<int> &,
//                         int, int);
//
// std::vector<int> max_card_sets(Domset&, int);
// bool p_cardinality(Domset&, int);
// int single_subsets(Domset&);
// bool included_sets(Domset&);
// int largest_logn(int);
// int largest_harmonic_num(int);
//
//
// //for TESTING
// bool is_cover(Domset &, Map<Set*>*, OPT_2_1_Sets&);
// bool is_domset(Graph*, Set*);
// void print_reduction(Domset&);
// void print_map(Map<Set*>*);
// void print_all(Domset &, Map<Set*>*, Map<Set*>*, Map<Set*>*);

#endif
