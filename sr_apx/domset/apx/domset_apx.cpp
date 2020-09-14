
#include <cstdio>
#include <algorithm>  //max element

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


void mod_exp_c_apx(Graph* graph, Set* domset, Set* removedC, Set* removedS, int q, int ITER) {
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
    printf("iter=%d\n", ITER);
    if(ITER==10) return;
    
    
    int p=largest_harmonic_num(q);
    printf("\n\nq=%d, p=%d\n", q, p);
    
    while(removedC->size() < graph->size()) {  //TODO when a vert is added to domset, nbs must be added to C.
        printf("domset size=%d, removedC size=%d, removedS size=%d\n", domset->size(), removedC->size(), removedS->size()); 
        
        printf("\ndomset: ");
        for(auto it=domset->begin(); it!=domset->end(); it++) printf(" %d,\n", *it);
        
        printf("\nremovedC: ");
        for(auto it=removedC->begin(); it!=removedC->end(); it++) printf(" %d,\n", *it);
        
        printf("\nremovedS: ");
        for(auto it=removedS->begin(); it!=removedS->end(); it++) printf(" %d,\n", *it);
        
        //1. IF there exists an item of C that belongs to a single subset S ∈ S, THEN add S to the solution
        single_subset(graph, domset, removedC, removedS);
        
        
        //2. IF there exist two sets S, R in S such that S is included into R, THEN remove S without branching;
        included_sets(graph, removedC, removedS);
        
        
        /* 3. IF all the residual subsets have cardinality at most p, 
         *    THEN run the algorithm by [Duh & Furer '97] in order to compute a q-approximation of the 
         *    optimal solution in the surviving instance
         */  
        bool run_apx = p_cardinality(graph, removedC, removedS, p);
        run_apx=false; 
        
        if(run_apx) {
            domset->add_all(Hk_minus_half_apx(graph, removedC, removedS, p)); //NOTE ?
        }
        
        /* 4. determine q sets S1, . . . , Sq from S such that (union i<=q Si) has maximum cardinality and perform 
        *    the the following branching: 
        * 
        *      a. either add every Si to the solution (and remove (union i<=q Si) from C), and remove the sets,
        *      b. or remove all of them.
        */ 
        
        std::vector<int> maxcard = max_card_sets(graph, removedC, removedS, q);
        printf("max card size=%d\n", maxcard.size());
        for(auto it=maxcard.begin(); it!=maxcard.end(); it++) printf("  maxcard verts=%d\n", *it);
        
        ITER=ITER+1;
        mod_exp_c_apx(graph, addto_set(domset, maxcard, false), addto_set(removedC, maxcard, false), 
                      addto_set(removedS, maxcard, false), q, ITER); 
        //mod_exp_c_apx(graph, domset, removedC, addto_set(removedS, maxcard), q);
        
     }
}


Set* Hk_minus_half_apx(Graph* graph, Set* removedC, Set* removedS, int k) {  
    /*
     * Approximation algorithm (originally for set cover) from "Approximationof k-SetCover by Semi-Local Optimization"
     * [Duh & Furer 1997].
     * 
     * This version of the function is used as a subroutine in mod_exp_c_apx();
     */
    int l;
    if(k>=5) l=5;
    else if(k==4) l=4;
    //else if(k<4) //TODO  --instance of 3-set cover
    
    printf("Hk alg!, k=%d, l=%d\n", k, l);
    Set* domset = new Set();
    Set* removedC_temp = new Set();
    Set* removedS_temp = new Set();
    
    //Greedy phase-greedily choose a maximal collection of j-sets / vertices of deg. j
    maximal_jsets(graph, domset, removedC, removedC_temp, removedS, removedS_temp, k, l);
    
    /* Restricted Phase: choose maximal collection of j sets w. restriction that the 
     * choice of these j-sets wont increase the number of 1-sets in the chosen solution,
     * 
     * 1-sets are vertices w. no neighbors (in this instance).
     */
    restricted_phase(graph, domset, removedC, removedC_temp, removedS, removedS_temp, k, l);
    
    
    /* Semi-local Improvement phase for 3-set cover: run semi-local optimization on 
     * the still uncovered elements of U.
     */
    semi_local_opt(graph, domset, removedC, removedC_temp, removedS, removedS_temp);
    
    
    //Add all *_temp to removedC/removedS
    removedC->add_all(removedC_temp);
    removedS->add_all(removedS_temp);
    
    return domset;
}




//-----helpers
void semi_local_opt(Graph* graph, Set* domset, Set* removedC, Set* removedC_temp, Set* removedS, Set* removedS_temp) {
    /* Semi-local Improvement phase for 3-set cover: run semi-local optimization on 
     * the still uncovered elements of U.
     * 
     * NOTE this step is the reason we keep a separate removedC_temp/removedS_temp (or it should be).
     */
    //TODO
}


void restricted_phase(Graph* graph, Set* domset, Set* removedC, Set* removedC_temp, Set* removedS, Set* removedS_temp, int k, int l) {
    for(int j=l; j>=4; j--) {
        for(auto it=graph->begin(); it!=graph->end(); it++) {
            int v=*it;
            
            //checks if its in the current instance && uncovered by this step currently.
            if(!removedS->contains(v) && !removedS_temp->contains(v)) {           
                int deg_v=vertex_degree(graph, removedC, removedC_temp, v);
                if(deg_v == j) {
                    if(!does_set_increase_onesets(graph, removedC, removedC_temp, removedS, removedS_temp)) {
                        domset->insert(v);
                        removedS_temp->insert(v);  //already used this set in solution. 
                        addto_set(removedC_temp, graph->neighbors(v), removedC); //neighbors are now covered
                    }
                }
            }
        }
    }
}

int number_one_sets(Graph* graph, Set* removedC, Set* removedC_temp, Set* removedS_temp) {
    /* Finds the number of one sets needed to construct the rest of the solution,
     * given the partial solution up to this point.
     */
    //TODO redo this
    int current_num_onesets=0;
    for(auto it=removedS_temp->begin(); it!=removedS_temp->end(); it++) {
        int deg = vertex_degree(graph, removedC, removedC_temp, *it);
        if(deg == 0) current_num_onesets++;
    }
    return current_num_onesets;
}


bool does_set_increase_onesets(Graph* graph, Set* removedC, Set* removedC_temp, Set* removedS, Set* removedS_temp) {
    //want this to be false.
    bool increases = false;
    int current_num_onesets = number_one_sets(graph, removedC, removedC_temp, removedS_temp);
    
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int v=*it;
        //TODO
    }
    
    return increases;
}


void maximal_jsets(Graph* graph, Set* domset, Set* removedC, Set* removedC_temp, Set* removedS, Set* removedS_temp, int k, int l) {
    //step 1 in hk apx set covering algorithm
    for(int j=k; j>=l+1; j--) {
        //first we must get the degree of the vertex
        for(auto it=graph->begin(); it!=graph->end(); it++) {
            int v=*it;
            
            //checks if its in the current instance && uncovered by this step currently.
            if(!removedS->contains(v) && !removedS_temp->contains(v)) {           
                int deg_v=vertex_degree(graph, removedC, removedC_temp, v);
                if(deg_v == j) {
                    domset->insert(v);
                    removedS_temp->insert(v);  //already used this set in solution. 
                    addto_set(removedC_temp, graph->neighbors(v), removedC); //neighbors are now covered
                }
            }
        }
    }
}

int vertex_degree(Graph* graph, Set* removedC, Set* removedC_temp, int vert) {
    int degree=0;
    for(auto it=graph->neighbors(vert)->begin(); it!=graph->neighbors(vert)->end(); it++) {
        if(!removedC->contains(*it) && !removedC_temp->contains(*it)) degree++;
    }
    return degree;
}


int vertex_degree(Graph* graph, Set* removedC, int vert) {
    int degree=0;
    for(auto it=graph->neighbors(vert)->begin(); it!=graph->neighbors(vert)->end(); it++) {
        if(!removedC->contains(*it)) degree++;
    }
    return degree;
}


void addto_set(Set* addingto, Set* adding, Set* removedC) {
    //adds elements in adding\removedC to addingto set.
    for(auto it=adding->begin(); it!=adding->end(); it++) {
        if(!removedC->contains(*it)) addingto->insert(*it);
    }
}


Set* addto_set(Set* set, std::vector<int> vertices, bool add_nbs) {  //TODO add neigbors to set also.
    for(auto it=vertices.begin(); it!=vertices.end(); it++) set->insert(*it);
    return set;
}


Set* union_Ssets(Graph* graph, Set* removedC, Set* removedS, std::vector<int> tmp) {
    Set* un = new Set();
    for(int i=0; i<tmp.size(); i++) {
        Set* nbs_u = graph->neighbors(tmp[i]);
        for(auto it=nbs_u->begin(); it!=nbs_u->end(); it++) {
            int v=*it;
            if(!removedC->contains(v) && !removedS->contains(v)) un->insert(v);
        }
        if(!removedC->contains(tmp[i]) && !removedS->contains(tmp[i])) un->insert(tmp[i]);
    }
    return un;
}


void vertex_combination(Graph* graph, Set* removedC, Set* removedS, 
                        std::vector<int> &ans, 
                        std::vector<int> &un_sizes,
                        std::vector<int> &vertices, 
                        std::vector<int> &tmp, 
                        int left, int q) {
    //tries all combinations of the sets, and finds set w/ max size of their union. O( q*n + q*delta ) time?    
    if(q==0) { //base case --q*delta time
        Set* un = union_Ssets(graph, removedC, removedS, tmp);
        
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
    
    for(int i=left; i<vertices.size(); i++) {
        tmp.push_back(vertices[i]);
        
        vertex_combination(graph, removedC, removedS, ans, un_sizes, vertices, tmp, i+1, q-1);
        tmp.pop_back();
    }
}


std::vector<int> max_card_sets(Graph* graph, Set* removedC, Set* removedS, int q) {
    /* Returns a set containing q vertices corresponding to 
     * 
     * Finds a set of vertices, not in removedC set, whose union of neighbors plus themselves is maximum.
     */    
    std::vector<int> vertices_left;
    for(auto it=graph->begin(); it!=graph->end(); it++) {  //n 
        int v=*it;
        if(!removedC->contains(v) && !removedS->contains(v)) {
            vertices_left.push_back(v);
        }
    }
    
    std::vector<int> max_card_comb;
    std::vector<int> un_sizes;
    std::vector<int> tmp;
    vertex_combination(graph, removedC, removedS, max_card_comb, un_sizes, vertices_left, tmp, 0, q);
    
    return max_card_comb;
}


bool p_cardinality(Graph* graph, Set* removedC, Set* removedS, int p) { 
    //Returns true IF all the residual subsets have cardinality at most p--degree plus themselves
    bool run_apx=true;
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int v=*it;
        if(!removedS->contains(v) && !removedC->contains(v)) {
            int sizeof_S = vertex_degree(graph, removedC, v) + 1;
            
            if(sizeof_S > p) return false;
        }
    }
    
    return run_apx;  //true
}


void single_subset(Graph* graph, Set* domset, Set* removedC, Set* removedS) {
    /*1. IF there exists an item of C=V\removedC that belongs to a single subset S ∈ S\removedS, 
     * THEN add S to the solution
     * 
     * For dom set version, a vertex belongs to only one subset if it is a single vertex w/ no neighbors.
     */
    printf("finding single subsets\n");
    
    for(auto it=graph->begin(); it!=graph->end(); it++) { //m*n
        int u=*it;
        
        if(!removedC->contains(u) && !removedS->contains(u)) { 
            int deg_u = vertex_degree(graph, removedC, u);
            
            if(deg_u==0) {
                domset->insert(u);  //add to soln
                removedC->insert(u); //now removed NOTE double check that we do this
            }
        }
    }
}


void included_sets(Graph* graph, Set* removedC, Set* removedS) {
    /* 2. IF there exist two sets S, R in S such that S is included into R, THEN remove S without branching;
     * 
     * For domset, check each vertex,
     */
    printf("fining included sets\n");
    for(auto it=graph->begin(); it!=graph->end(); it++) { 
        int u=*it;
    
        if(!removedS->contains(u) && !removedC->contains(u)) {
            Set* nbs_u = graph->neighbors(u);  //check nbs_u+u subset of nbs_v+v
            for(auto itt=nbs_u->begin(); itt!=nbs_u->end(); itt++) { 
                bool is_subset=true;
                int v=*itt;
                
                if(!removedS->contains(v) && !removedC->contains(v)) {
                    Set* nbs_v = graph->neighbors(v);
                    for(auto inb=nbs_u->begin(); inb!=nbs_u->end(); inb++) { 
                        int nb_u=*inb;
                        if(!nbs_v->contains(nb_u) && v!=nb_u && !removedC->contains(nb_u)) {
                            is_subset=false; 
                        }
                    }
                    
                    if(is_subset) {
                        removedS->insert(u);  //only insert the vertex, so we know not to look at its neighbors in future.
                    }
                }
            }
        }
    }
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


//For TESTING
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
