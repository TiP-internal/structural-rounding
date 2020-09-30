
#include <cstdio>
#include <algorithm>  //max element
#include <math.h>       /* exp */

#include "domset_apx.hpp"


Set* logn_apx(Graph* graph) {
    /*
     * Chooses vertex of max degree and adds to dominating set.
     * Adds all of the neighbors and itself to the visited set.
     * Repeats while some vertices have not been visited yet.
     */
    Set* domset = new Set();
    Set* visited = new Set();

    while(visited->size() < graph->size()) {
        int max_deg_v = max_deg_vertex(graph, visited);

        domset->insert(max_deg_v);
        visited->insert(max_deg_v);
        for(auto it=graph->neighbors(max_deg_v)->begin(); it!=graph->neighbors(max_deg_v)->end(); it++) {
            visited->insert(*it);
        }
    }
    delete visited;
    return domset;
}

int logn_apx(Domset &domset_inst) {
    //NOTE this is for is we want to make the mod exp algorithm
    //more dom set specific

    Set* domset = new Set();

    return 0;
}


Domset reduction(Graph* graph) {
    /* Reduces a dominating set instance to a set cover instance.
     */
    Domset ds_reduct;
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int v=*it;
        ds_reduct.base->insert(v);

        Set* nbs_set = new Set();
        nbs_set->insert(v);
        for(auto itt=graph->neighbors(v)->begin(); itt!=graph->neighbors(v)->end(); itt++) {
            int u=*itt;
            nbs_set->insert(u);
        }
        ds_reduct.collection->insert(v, nbs_set); //key, value?
    }
    return ds_reduct;
}


Set* mod_exponential_domset(Domset &domset_inst, int q) {
    int size = mod_exp_c_apx(domset_inst, q);
    delete domset_inst.base, domset_inst.collection;
    return domset_inst.domset;
}


int mod_exp_c_apx(Domset domset_inst, int q) {
    /*
     * Approximately prunes search tree.
     * Gives q-approximation in O∗((α1 * α2)^n ) time.
     *
     * m=number of sets
     * Proposition 2. For any integer q >= 1, Algorithm SC1 computes a q-approximation of min set cover in O∗(2^m/q )
     *
     * Proposition 6. Assume there exists an r-approximation algorithm A for min set cover (r >= 1)
     * with running time O∗(α^n_1 * α^m_2 ). Then, there exists an r-approximation algorithm for min dominating
     * set with running time O∗((α1 * α2)^n ).
     *
     * The reduction from set cover to dominating set:
     * "Let G(V, E) be an instance of min dominating set. We construct an instance I(S, C) of min set cover
     * as follows: C = V, S = {S_v = {v} ∪ Γ(v), v ∈ V}, where Γ(v) is the set of neighbors of vertex v (|S| = |V|).
     * Consider now a cover S' = {S_v1, . . . , S_vk } of C. Obviously, the set {v_1, . . . , v_k} is a
     * dominating set of G, since set S_vi (resp., vertex vi) covers (resp., dominates) elements
     * corresponding to vertex vi itself and to its adjacent vertices."
     */
    if(domset_inst.base->size()==0) return 0;

    int p=largest_harmonic_num(q);

    /* 1. IF there exists an item of C that belongs to a single subset S ∈ S,
     * THEN add S to the solution
     */
    int val = single_subsets(domset_inst);
    if(val>0) return val + mod_exp_c_apx(domset_inst, q);

    /* 2. IF there exist two sets S, R in S such that S is included into R,
     * THEN remove S without branching;
     */
    bool included = included_sets(domset_inst);

    if(included) return mod_exp_c_apx(domset_inst, q);

    /* 3. IF all the residual subsets have cardinality at most p,
     *    THEN run the algorithm by [Duh & Furer '97] in order to compute a
     *    q-approximation of the
     *    optimal solution in the surviving instance
     */
    bool run_apx = p_cardinality(domset_inst, p);
    if(run_apx) {
        int apx = Hk_minus_half_apx(domset_inst, p);
        return apx + mod_exp_c_apx(domset_inst, q);
    }

    /* 4. determine q sets S1, . . . , Sq from S such that (union i<=q Si)
        * has maximum cardinality and perform
    *    the the following branching:
    *
    *      a. either add every Si to the solution (and remove
    *         (union i<=q Si) from C), and remove the sets,
    *      b. or remove all of them.
    */
    std::vector<int> maxcard;
    if(q <= domset_inst.collection->size()) {
        maxcard = max_card_sets(domset_inst, q);
    } else {
        for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
            maxcard.push_back(*it);
        }
    }

    int count = how_many(domset_inst, maxcard);
    return std::min(mod_exp_c_apx(remove_sets_from_S(domset_inst, maxcard), q),
             (count+mod_exp_c_apx(remove_from_SandC_addsoln(domset_inst, maxcard), q)));
}


int Hk_minus_half_apx(Domset &domset_inst, int k) {
    /*
     * Approximation algorithm (originally for set cover) from
     * "Approximationof k-SetCover by Semi-Local Optimization"
     * [Duh & Furer 1997].
     *
     * This version of the function is used as a subroutine in mod_exp_c_apx();
     */
    int l;
    if(k>=5) l=5;
    else if(k==4) l=4;
    else if(k<4) l=3;   //--instance of 3-set cover, only run semi_local_opt

    int num_large_sets_added = 0;  //number of >3 sized sets added to soln

    //Greedy phase-greedily choose a maximal collection of j-sets (sets of size >= 6)
    num_large_sets_added = maximal_jsets(domset_inst, k, l);

    /* Restricted Phase: choose maximal collection of j sets w. restriction that the
     * choice of these j-sets wont increase the number of 1-sets in the chosen solution,
     * Finds 5 and/or 4 sets only in this phase.
     */
    num_large_sets_added += restricted_phase(domset_inst, k, l);

    /* Semi-local Improvement phase for 3-set cover: run semi-local optimization on
     * the still uncovered elements of U.
     */
    num_large_sets_added+=semi_local_opt(domset_inst, num_large_sets_added);

    return num_large_sets_added;
}



//-----helpers
int semi_local_opt(Domset &domset_inst, int num_large_sets_added) {
    /* Semi-local Improvement phase for 3-set cover: run semi-local optimization on
     * the still uncovered elements of U.
     *
     * Semi-local optimization can easily be explained in the context of 3-Set Cover.
     * In a pure local improvement step, a current approximate solution is improved by
     * replacing a constant number of sets in the cover with a (hopefully smaller) set
     * of other admissible small sets to obtain a new cover. For semi-local improvements,
     * we observe that once the sets of size 3 have been selected for a partial cover,
     * then opt_2_1_set_cover finds the remainder of 2/1 sets optimally (global changes
     * possible).
     *
     * Hence. a semi-local (s, t)-improvement step for 3-Set Cover consists of the insertion
     * of up to s 3-sets, and the deletion of up to t 3-sets from a current cover. In
     * addition, an arbitrary number of 2-sets and 1-sets are optimally replaced.
     *
     * This means that we only count the number of 3-sets involved instead of all the
     * small subsets.
     *
     * This scheme thus allows globally changes of 2-sets and 1-sets.
     *
     * The objective function used to measure improvements is not necessary equal to the
     * given objective function. (i.e., the number of sets in the cover.) We will order
     * the quality of a solution lexicographically, first by the number of sets in the cover,
     * and second by the number of 1-sets, i.e., smaller covers are better, and among covers
     * of a given size, we prefer those with fewer 1-sets.
     *
     *
     * (2,1)-improvements specifically.
     *
     * We are given a partial solution up to this point, covered by 4-sets, 5-sets, ... k-sets.
     * In this step we need to use (2,1)-improvements to find the best 3-sets,
     * Once the 3-sets are chosen, find the best 2/1-sets optimally.
     *
     * Should try to find local improvements until solution stops improving.
     */
    int num_sets_added=0;
    Map<Set*>* three_sets_added = new Map<Set*>();

    //Step 1: Greedily select a maximal set of 3-sets
    greedy_3set_select(domset_inst, three_sets_added);

    //Step 3: local (2,1)-improvements
    num_sets_added = st_improvements(domset_inst, three_sets_added, num_large_sets_added);
    return num_sets_added;
}


OPT_2_1_Sets opt_2_1_set_cover(Domset &domset_inst, Map<Set*>* three_sets_added) {
    /* the task of covering the remainder by 2-sets and 1-sets can be done optimally
     * in polynomial time by computing a maximum matching [6]. The vertices are the
     * elements of U still to cover and the edges are the admissible 2-sets. The
     * 2-sets corresponding to the maximum matching edges and the 1-sets corresponding
     * to the vertices not covered by the maximum matching edges form a optimum
     * covering of the 2-Set Cover,
     */
    OPT_2_1_Sets sets;

    //first create auxiliary graph-edge if there is a 2-set between 2 uncovered vertices
    Graph* aux_graph = new Graph();

    std::vector<int> keys_toremove;
    for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
        int key = *it;
        Set* curr_set = domset_inst.collection->at(key);

        //NOTE finds set size excluding the keys of the three sets
        int set_size=0;
        int u=-99, v=-99;
        for(auto itt=curr_set->begin(); itt!=curr_set->end(); itt++) {
            int val=*itt;
            bool good=true;
            for(auto itset=three_sets_added->begin(); itset!=three_sets_added->end(); itset++) {
               int k = *itset;
               Set* s = three_sets_added->at(k);

               if(s->contains(val)) good=false; //val is actually covered by a set in three sets
               else {
                    if(u==-99) u=val;
                    else if(v==-99) v=val;
               }
           }

           if(good) set_size++;
        }

        if(set_size == 2) {
            aux_graph->add_edge(u, v);
        }

        if(set_size == 1) {
            //If theres a 1-set, add to the solution.
            sets.one_sets_added->insert(key, curr_set);
            keys_toremove.push_back(key);
        }
    }

    for(auto it=aux_graph->begin(); it!=aux_graph->end(); it++) {
        int v=*it;
        for(auto itt=aux_graph->neighbors(v)->begin();
            itt!=aux_graph->neighbors(v)->end(); itt++) {
            int u=*itt;
        }
    }

    for(int i=0; i<keys_toremove.size(); i++) domset_inst.collection->remove(keys_toremove[i]);

    //now find maximum matching? -- this is supposed to be done optimally.
    max_matching(aux_graph, domset_inst, sets, three_sets_added);

    delete aux_graph;

    return sets;
}


std::vector<std::vector<int>> threeset_combs(Map<Set*>* three_sets_added) {
    std::vector<std::vector<int>> combos;

    //2 or 1 3-sets get replaced by 0 of 1 different 3-sets
    for(auto it=three_sets_added->begin(); it!=three_sets_added->end(); it++) {
        int key = *it;
        for(auto itt=three_sets_added->begin(); itt!=three_sets_added->end(); itt++) {
            int key_inner = *itt;
            std::vector<int> comb;
            comb.push_back(key);
            comb.push_back(key_inner);

            combos.push_back(comb);
        }
    }
    return combos;
}


int st_improvements(Domset &domset_inst, Map<Set*>* three_sets_added,
                    const int num_large_sets_added) {
    /* (s,t)-improvement:
     *
     * take current (full) apx soln,
     * improve it by replacing up to 2 sets w/ up to 1 such set.
     * (e.g. replace up to 2 3-sets w/ up to 1 3-set)
     *
     * --replace (0,1,or 2) 3-sets w/ (0, or 1) 3-sets.
     * --the solution is improved if below criteria met.
     * --if improved take the new solution, if not improved dont take new solution.
     *
     * an improvement occurs when the number of sets in the cover stays the same or
     * decreases, and/or the number of 1-sets decreases. --fewer 1-sets are better.
     *
     * cases (in which to accept the updated solution):
     *  a=|S| -> size of current solution
     *  b=|1-sets| -> number of one sets in the solution
     *
     *      a. a stays the same size, and b decreases
     *      b. a decreases, b stays the same, decreases, (and when b increases?)
     *
     * NOTE this function could be cleaned up a bit
     */
    OPT_2_1_Sets sets = opt_2_1_set_cover(domset_inst, three_sets_added);

    //---adds 1/2 sets back to collection after removal in opt_2_1_set_cover step
    for(auto it=sets.two_sets_added->begin(); it!=sets.two_sets_added->end(); it++) {
        int key=*it;
        Set* curr_set = sets.two_sets_added->at(key);
        if(!domset_inst.collection->contains(key)) {
            domset_inst.collection->insert(key, curr_set);
        }
    }

    for(auto it=sets.one_sets_added->begin(); it!=sets.one_sets_added->end(); it++) {
        int key=*it;
        Set* curr_set = sets.one_sets_added->at(key);
        if(!domset_inst.collection->contains(key)) {
            domset_inst.collection->insert(key, curr_set);
        }
    }
    //---

    int a_current, b_current, a_new, b_new;
    b_current = sets.one_sets_added->size();
    //current solution size
    a_current = num_large_sets_added+three_sets_added->size()+sets.two_sets_added->size()+b_current;
    std::vector<std::vector<int>> combos = threeset_combs(three_sets_added);

    //Step 1: (2, 1)-improvements
    bool still_improving = true;
    while(still_improving) {
        //printf("______________________________________________________improving\n");
        bool improved=false;
        for(int i=0; i<combos.size(); i++) { //for each combination of 2 sets to replace
            Set* firstset = three_sets_added->at(combos[i][0]);
            Set* secset = three_sets_added->at(combos[i][1]);

            three_sets_added->remove(combos[i][0]); //remove from chosen 3-sets
            three_sets_added->remove(combos[i][1]);
            domset_inst.collection->insert(combos[i][0], firstset); //add back to the collection
            domset_inst.collection->insert(combos[i][1], secset);

            //----------replace w/ 1 different 3-sets
            bool replaced=false;
            for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
                int key = *it;
                Set* replacement_set = domset_inst.collection->at(key);

                if(replacement_set->size() == 3 && key!=combos[i][0] && key!=combos[i][1]) {
                    replaced=true;
                    three_sets_added->insert(key, replacement_set);
                    domset_inst.collection->remove(key);

                    OPT_2_1_Sets sets_temp = opt_2_1_set_cover(domset_inst, three_sets_added);
                    b_new = sets_temp.one_sets_added->size();
                    a_new = num_large_sets_added+three_sets_added->size()+sets_temp.two_sets_added->size()+b_new;

                    for(auto it=sets_temp.two_sets_added->begin(); it!=sets_temp.two_sets_added->end(); it++) {
                        int key=*it;
                        Set* curr_set = sets_temp.two_sets_added->at(key);
                        if(!domset_inst.collection->contains(key)) {
                            domset_inst.collection->insert(key, curr_set);
                        }
                    }

                    for(auto it=sets_temp.one_sets_added->begin(); it!=sets_temp.one_sets_added->end(); it++) {
                        int key=*it;
                        Set* curr_set = sets_temp.one_sets_added->at(key);
                        if(!domset_inst.collection->contains(key)) {
                            domset_inst.collection->insert(key, curr_set);
                        }
                    }

                    //decide whether to keep new solution or not.
                    if((a_current == a_new && b_new < b_current)
                        || (a_new < a_current)) {
                        //accept new solution
                        delete sets.two_sets_added, sets.one_sets_added;
                        sets.two_sets_added = sets_temp.two_sets_added;
                        sets.one_sets_added = sets_temp.one_sets_added;
                        delete sets_temp.two_sets_added, sets_temp.one_sets_added;

                        improved=true;
                        break;
                    } else {
                        //remove the just added set
                        three_sets_added->remove(key);
                        domset_inst.collection->insert(key, replacement_set);

                        //add sets back to solution
                        three_sets_added->insert(combos[i][0], firstset);
                        three_sets_added->insert(combos[i][1], secset);
                        domset_inst.collection->remove(combos[i][0]);
                        domset_inst.collection->remove(combos[i][1]);

                        delete sets_temp.two_sets_added, sets_temp.one_sets_added;
                    }
                }
            }
            if(!replaced) {
                //the removed 3-sets never got replaced, so we are adding them back here.
                //add sets back to solution
                three_sets_added->insert(combos[i][0], firstset);
                three_sets_added->insert(combos[i][1], secset);
                domset_inst.collection->remove(combos[i][0]);
                domset_inst.collection->remove(combos[i][1]);
            }
        }

        if(!improved) still_improving=false; //end while loop
    }

    //Here update domset_inst
    int num_sets_added=0;
    for(auto it=three_sets_added->begin(); it!=three_sets_added->end(); it++) {
        int key=*it;
        Set* curr_set = three_sets_added->at(key);
        three_sets_added->remove(key);

        if(!domset_inst.domset->contains(key)) {
            domset_inst.domset->insert(key);
            num_sets_added++;
        }
        domset_inst.collection->insert(key, curr_set); //so its in the collection
        remove_from_collection(domset_inst, key);
    }

    for(auto it=sets.two_sets_added->begin(); it!=sets.two_sets_added->end(); it++) {
        int key=*it;
        Set* curr_set = sets.two_sets_added->at(key);
        sets.two_sets_added->remove(key);

        if(!domset_inst.domset->contains(key)) {
            domset_inst.domset->insert(key);
            num_sets_added++;
        }
        remove_from_collection(domset_inst, key);
    }

    for(auto it=sets.one_sets_added->begin(); it!=sets.one_sets_added->end(); it++) {
        int key=*it;
        Set* curr_set = sets.one_sets_added->at(key);
        sets.one_sets_added->remove(key);

        if(!domset_inst.domset->contains(key)) {
            domset_inst.domset->insert(key);
            num_sets_added++;
        }
        remove_from_collection(domset_inst, key);
    }

    delete three_sets_added, sets.two_sets_added, sets.one_sets_added;
    return num_sets_added;
}


void greedy_3set_select(Domset &domset_inst, Map<Set*>* three_sets_added) {
    /*
     * Greedily select a maximal set of 3-sets???
     */
    //maximal_set(domset_inst, 3);

    for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
        int key = *it;
        Set* curr_set = domset_inst.collection->at(key);
        int set_size = setsize(three_sets_added, curr_set);

        if(set_size == 3) {
            three_sets_added->insert(key, curr_set);
            domset_inst.collection->remove(key);                //remove set from collection
        }
    }
}


int restricted_phase(Domset &domset_inst, int k, int l) {
    /* Restricted Phase: choose maximal collection of j sets w. restriction that the
     * choice of these j-sets wont increase the number of 1-sets in the chosen solution,
     * Here, j is either 5 and/or 4.
     */
    int num_sets = 0;
    for(int j=l; j>=4; j--) {
        for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
            int key = *it;
            Set* curr_set = domset_inst.collection->at(key);

            if(curr_set->size() == j) {
                //does deleting key increase num onesets?
                bool incr_onesets = does_del_incr_onesets(domset_inst, key);
                if(!incr_onesets) { //if not go ahead
                    num_sets++;
                    domset_inst.domset->insert(key);                    //add to solution
                    remove_from_collection(domset_inst, key);           //remove value from all other sets
                }
            }
        }
    }
    return num_sets;
}


Domset remove_sets_from_S(Domset domset_inst, std::vector<int> maxcard_sets) {
    for(int i=0; i<maxcard_sets.size(); i++) {
        domset_inst.collection->remove(maxcard_sets[i]);
    }
    return domset_inst;
}


int how_many(Domset &domset_inst, std::vector<int> &maxcard_sets) {
    int count=0;
    for(int i=0; i<maxcard_sets.size(); i++) {
        if(!domset_inst.domset->contains(maxcard_sets[i])) count++;
    }
    return count;
}


Domset remove_from_SandC_addsoln(Domset domset_inst, std::vector<int> maxcard_sets) {
    for(int i=0; i<maxcard_sets.size(); i++) {
        //remove from C
        for(auto it=domset_inst.collection->at(maxcard_sets[i])->begin();
            it!=domset_inst.collection->at(maxcard_sets[i])->end(); it++) {
            domset_inst.base->remove(*it);
        }
        //add to solution
        domset_inst.domset->insert(maxcard_sets[i]);

        //remove from S
        domset_inst.collection->remove(maxcard_sets[i]);
    }
    return domset_inst;
}


int maximal_jsets(Domset &domset_inst, int k, int l) {
    /* step 1 in hk apx set covering algorithm
     * Greedy phase-greedily chooses a maximal collection of j-sets
     */
    int num_sets=0;
    for(int j=k; j>=l+1; j--) {
        num_sets += maximal_set(domset_inst, j);
    }
    return num_sets;
}


int maximal_set(Domset &domset_inst, int j) {
    //Finds sets of size j, removes them from colleciton, adds key to domset solution.
    int num_sets=0;
    for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
        int key = *it;
        Set* curr_set = domset_inst.collection->at(key);

        if(curr_set->size() == j) {
            num_sets++;
            domset_inst.domset->insert(key);                    //add to solution
            //NOTE remove value from all other sets
            remove_from_collection(domset_inst, key);
        }
    }
    return num_sets;
}


void vertex_combination(Domset &domset_inst,
                        std::vector<int> &ans,
                        std::vector<int> &un_sizes,
                        std::vector<int> &tmp,
                        int left, int q) {
    //tries all combinations of the sets, and finds set w/ max size of their union. O( q*n + q*delta ) time?
    if(q==0) { //base case --q*delta time
        Set* un = new Set();
        for(auto it=tmp.begin(); it!=tmp.end(); it++) {
            Set* adding = domset_inst.collection->at(*it);
            un->add_all(adding);
        }

        int maxun;
        if(un_sizes.size() > 0) {
            maxun=*max_element(std::begin(un_sizes), std::end(un_sizes));
        } else maxun=0;

        if(un->size() > maxun) {
            maxun = un->size();
            ans = tmp;
        }
        un_sizes.push_back(un->size());

        delete un;
        return;
    }

    int i=0;
    auto it=domset_inst.collection->begin();
    for(int i=0; i<left; i++) it++;  //sets it=left

    i=left;
    for(it; it!=domset_inst.collection->end(); it++) {
        tmp.push_back(*it);
        vertex_combination(domset_inst, ans, un_sizes, tmp, i+1, q-1);
        tmp.pop_back();

        i++;
    }
}


std::vector<int> max_card_sets(Domset &domset_inst, int q) {
    /* determine q sets S1, . . . , Sq from S such that (union i<=q Si) has maximum cardinality and perform
     *
     * Returns a set containing q vertices corresponding to
     *
     * Finds a set of vertices, not in removedC set, whose union of neighbors plus themselves is maximum.
     */
    std::vector<int> ans;
    std::vector<int> un_sizes;
    std::vector<int> tmp;

    vertex_combination(domset_inst, ans, un_sizes, tmp, 0, q);

    return ans;
}


int setsize(Map<Set*>* three_sets_added, Set* curr_set) {
    //---gets set size
    int set_size=0;
    for(auto itt=curr_set->begin(); itt!=curr_set->end(); itt++) {
        int val=*itt;
        bool good=true;
        for(auto itset=three_sets_added->begin(); itset!=three_sets_added->end(); itset++) {
            int k = *itset;
            Set* s = three_sets_added->at(k);
            if(s->contains(val)) good=false; //val is actually covered by a set in three sets
        }

        if(good) set_size++;
    }
    return set_size;
}


void max_matching(Graph* aux_graph, Domset &domset_inst,
                  OPT_2_1_Sets &sets, Map<Set*>* three_sets_added) {
    /* Needs to be implementation from:
     *
     * NOTE: But I dont have access to paper, so a greedy routine is used
     * temporarily
     */
    for(auto it=aux_graph->begin(); it!=aux_graph->end(); it++) {
        int v=*it;
        Set* nbs_v = aux_graph->neighbors(v);
        if(nbs_v->size()==0) { //add v to soln?
            sets.one_sets_added->insert(v, domset_inst.collection->at(v));
            domset_inst.collection->remove(v);

        } else {
            //put 2-set consisting of {v, nb} in solution
            int nb = *nbs_v->begin(); //select a neighbor
            for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
                int key = *it;
                Set* curr_set = domset_inst.collection->at(key);
                int set_size = setsize(three_sets_added, curr_set);
                if(set_size==2 && curr_set->contains(v) && curr_set->contains(v)) {
                    sets.two_sets_added->insert(key, curr_set);
                    domset_inst.collection->remove(key);
                }
            }
            //remove v+nbs from graph
            aux_graph->remove_vertex(nb);
            aux_graph->remove_vertex(v);
        }
    }
}


bool p_cardinality(Domset &domset_inst, int p) {
    //Returns true IF all the residual subsets have cardinality at most p--degree plus themselves
    bool run_apx=true;
    for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
        int key = *it;
        Set* curr_set = domset_inst.collection->at(key);

        if(curr_set->size() > p) run_apx=false;
    }

    return run_apx;  //true
}


int single_subsets(Domset &domset_inst) {
    /*1. IF there exists an item of C that belongs to a single subset S ∈ S,
     * THEN add S to the solution
     *
     * Find sets of size 1 and add to solution
     *
     * For dom set version, a vertex belongs to only one subset if it is a single vertex w/ no neighbors.
     */
    int sets_added=0;
    for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
        int key = *it;
        Set* curr_set = domset_inst.collection->at(key);

        if(curr_set->size()==1) {
            domset_inst.collection->remove(key);                //NOTE remove or erase here??
            domset_inst.domset->insert(key);                    //adding to solution
            domset_inst.base->erase(key);
            sets_added++;
        }
    }
    return sets_added;
}


bool included_sets(Domset &domset_inst) {
    /* 2. IF there exist two sets S, R in S such that S is included into R, THEN remove S without branching;
     *
     * For domset, check each vertex,
     */
    bool applied=false;
    for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
        int key = *it;
        Set* curr_set = domset_inst.collection->at(key);

        for(auto itt=domset_inst.collection->begin(); itt!=domset_inst.collection->end(); itt++) {
            int key_check = *itt;
            if(key!=key_check) {
                Set* check_set = domset_inst.collection->at(key_check);

                //is check_set subset of curr_set??
                bool subset = is_subset(check_set, curr_set);
                if(subset) {
                    domset_inst.collection->remove(key_check);
                    applied=true;
                }
            }
        }
    }
    return applied;
}


int vertex_degree(Graph* graph, Set* removedC, int vert) {
    int degree= removedC->contains(vert) ? 0 : 1;
    for(auto it=graph->neighbors(vert)->begin(); it!=graph->neighbors(vert)->end(); it++) {
        if(!removedC->contains(*it)) degree++;
    }
    return degree;
}


int max_deg_vertex(Graph* graph, Set* visited) {
    //finds the maximum degree not in the visited set
    int max_deg_value = 0;
    int max_deg_vert=0;

    for (auto iu = graph->begin(); iu!=graph->end(); iu++) {
        int u = *iu;
        if(!visited->contains(u)) {
            int u_deg = vertex_degree(graph, visited, u);

            if ( u_deg > max_deg_value) {
                max_deg_value = u_deg;
                max_deg_vert = u;
            }
        }
    }
    return max_deg_vert;
}


int largest_logn(int q) {
    /* finds ln(?)+1 <= q
     */
    if(q==0) return 0;

    return exp (q-1);  //NOTE floor maybe?
}


int largest_harmonic_num(int q) {
    /*
     * fix q ∈ N∗ and compute the largest integer p such that H(p) − 1/2 <= q,
     * where H is the harmonic number sequence
     */
    if(q==0) return 0;

    float p=1;
    float harm_num=0;
    while((harm_num-0.5) <= q) {
        harm_num+=1.0/p;
        p++;
    }
    return p-2;
}


bool is_subset(Set* A, Set* B) {
    //is A a subset of B?
    bool subset=true;
    for(auto it=A->begin(); it!=A->end(); it++) {
        if(!B->contains(*it)) subset=false;
    }
    return subset;
}


void remove_from_collection(Domset &domset_inst, int v) {
    //removes the vertex v from each set in the collection of sets.
    if(domset_inst.collection->contains(v)) {
        Set* removingset = domset_inst.collection->at(v);
        domset_inst.collection->remove(v);
        for(auto itt=removingset->begin(); itt!=removingset->end(); itt++) {
            int removing_val = *itt;
            domset_inst.base->remove(removing_val);

            for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
                int key = *it;
                Set* curr_set = domset_inst.collection->at(key);

                if(curr_set->contains(removing_val)) curr_set->remove(removing_val);  //NOTE remove or erase??
                if(curr_set->size()==0) domset_inst.collection->remove(key);
            }
        }
        delete removingset;
    }
}


bool does_del_incr_onesets(Domset &domset_inst, int key) {
    /* Does deleting all values in the key set increase single sets in collection?
     *
     * EDIT number of single sets necessary to form a cover.
     */

    bool increases = false;
    Set* removingset = domset_inst.collection->at(key);

    for(auto it=domset_inst.collection->begin(); it!=domset_inst.collection->end(); it++) {
        int keyy = *it;
        int num_removed=0;
        int value_left=-99;
        Set* curr_set = domset_inst.collection->at(keyy);

        for(auto itt=removingset->begin(); itt!=removingset->end(); itt++) {
            int removing_val = *itt;
            if(curr_set->contains(removing_val)) num_removed++;
            else value_left=removing_val;
        }

        if(curr_set->size()-num_removed == 1) {
            //check if the value in the potential one set is in any other sets. it it is then the 4/5 set is fine
            bool contained=false;
            for(auto ittt=domset_inst.collection->begin(); ittt!=domset_inst.collection->end(); ittt++) {
                int keyyy = *ittt;
                Set* s = domset_inst.collection->at(keyyy);
                if(s->contains(value_left)) contained=true;
            }

            if(!contained) increases=true;
        }
    }

    return increases;
}



//For TESTING
bool is_cover(Domset &domset_inst, Map<Set*>* three_sets_added, OPT_2_1_Sets &sets) {
    bool is_cover=true;

    for(auto it=domset_inst.base->begin(); it!=domset_inst.base->end(); it++) {
        int value=*it;
        bool covered=false;

        for(auto itt=three_sets_added->begin(); itt!=three_sets_added->end(); itt++) {
            int key = *itt;
            Set* curr_set = three_sets_added->at(key);
            if(curr_set->contains(value)) covered=true;
        }

        if(!covered) {
            for(auto itt=sets.two_sets_added->begin(); itt!=sets.two_sets_added->end(); itt++) {
                int key = *itt;
                Set* curr_set = three_sets_added->at(key);
                if(curr_set->contains(value)) covered=true;
            }
        }

        if(!covered) {
            for(auto itt=sets.one_sets_added->begin(); itt!=sets.one_sets_added->end(); itt++) {
                int key = *itt;
                Set* curr_set = three_sets_added->at(key);
                if(curr_set->contains(value)) covered=true;
            }
        }

        if(!covered) is_cover=false;
    }

    return is_cover;
}

bool is_domset(Graph* graph, Set* domset) {
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int v=*it;
        bool adjacent = false;

        for(auto ids=domset->begin(); ids!=domset->end(); ids++) {
            int u=*ids;
            if(graph->adjacent(v, u)) {
                adjacent=true;
            }
        }
        if(!adjacent && !domset->contains(v)) return false;
    }
    return true;
}


void print_reduction(Domset &ds) {
    printf("_________________________Printing Reduction\n");

    printf("Base Set: ");
    for(auto it=ds.base->begin(); it!=ds.base->end(); it++) printf(" %d,", *it);

    printf("\nCollection: \n");
    for(auto it=ds.collection->begin(); it!=ds.collection->end(); it++) {
        int key = *it;
        Set* col = ds.collection->at(key);
        printf(" \nkey: %d  |  vals: ", key);

        for(auto it=col->begin(); it!=col->end(); it++) printf(" %d,", *it);
    }

    printf("\n-----Dom Set: ");
    for(auto it=ds.domset->begin(); it!=ds.domset->end(); it++) printf(" %d,", *it);
    printf("\n______________________________\n\n");
}


void print_map(Map<Set*>* map) {
    printf("_________________________Printing Map set\n");
    for(auto it= map->begin(); it!=map->end(); it++) {
        int key = *it;
        Set* col = map->at(key);
        printf(" \nkey: %d  |  vals: ", key);

        for(auto it=col->begin(); it!=col->end(); it++) printf(" %d,", *it);
    }
    printf("\n______________________________\n\n");
}


void print_all(Domset &ds, Map<Set*>* three_map, Map<Set*>* two_map, Map<Set*>* one_map) {
    printf("\n____________________________________________________\n");

    printf("Base Set: ");
    for(auto it=ds.base->begin(); it!=ds.base->end(); it++) printf(" %d,", *it);

    printf("\nCollection: \n");
    for(auto it=ds.collection->begin(); it!=ds.collection->end(); it++) {
        int key = *it;
        Set* col = ds.collection->at(key);
        printf(" key: %d  |  vals: ", key);

        for(auto it=col->begin(); it!=col->end(); it++) printf(" %d,", *it);
        printf("\n");
    }

    printf("-----Dom Set: ");
    for(auto it=ds.domset->begin(); it!=ds.domset->end(); it++) printf(" %d,", *it);


    printf("\nThree Sets:\n");
    for(auto it= three_map->begin(); it!=three_map->end(); it++) {
        int key = *it;
        Set* col = three_map->at(key);
        printf(" key: %d  |  vals: ", key);

        for(auto it=col->begin(); it!=col->end(); it++) printf(" %d,", *it);
        printf("\n");
    }


    printf("Two Sets:\n");
    for(auto it= two_map->begin(); it!=two_map->end(); it++) {
        int key = *it;
        Set* col = two_map->at(key);
        printf(" key: %d  |  vals: ", key);

        for(auto it=col->begin(); it!=col->end(); it++) printf(" %d,", *it);
        printf("\n");
    }


    printf("One Sets:\n");
    for(auto it= one_map->begin(); it!=one_map->end(); it++) {
        int key = *it;
        Set* col = one_map->at(key);
        printf(" key: %d  |  vals: ", key);

        for(auto it=col->begin(); it!=col->end(); it++) printf(" %d,", *it);
        printf("\n");
    }

    printf("____________________________________________________\n\n");
}
