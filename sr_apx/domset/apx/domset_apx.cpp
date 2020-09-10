
#include <cstdio>

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


Set* mod_exp_c_apx(Graph* graph, int q) {
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
    Set* domset=new Set();
    Set* removedS = new Set(); //in set cover terms, these are the subsets from S which are removed
    Set* removedC = new Set(); //these are the removed vertices from the C (aka V)
    
    int p=largest_harmonic_num(q);
    while(removedC->size() < graph->size()) {
        
        //1. IF there exists an item of C that belongs to a single subset S ∈ S, THEN add S to the solution
        single_subset(graph, domset, removedC, removedS);
        
        
        //2. IF there exist two sets S, R in S such that S is included into R, THEN remove S without branching;
        included_sets(graph, removedC, removedS);
        
        
        /* 3. IF all the residual subsets have cardinality at most p, 
         *    THEN run the algorithm by [Duh & Furer '97] in order to compute a q-approximation of the 
         *    optimal solution in the surviving instance
         */  
        bool run_apx = p_cardinality(graph, removedC, removedS, p);
        
        if(run_apx) {
            domset->add_all(Hk_minus_half_apx(graph, removedC, removedS)); //NOTE ?
        } else {  //NOTE Im not sure about this if else stmnt.
        
            /* 4. determine q sets S1, . . . , Sq from S such that (union i<=q Si) has maximum cardinality and perform 
            *    the the following branching: 
            * 
            *      a. either add every Si to the solution (and remove (union i<=q Si) from C), and remove the sets,
            *      b. or remove all of them.
            */ 
            
            //TODO
        }
        
    }
    
    delete removedC, removedS;
    return domset;
}


Set* Hk_minus_half_apx(Graph* graph, Set* removedC, Set* removedS) {  
    /*
     * Approximation algorithm (originally for set cover) from "Approximationof k-SetCover by Semi-Local Optimization"
     * [Duh & Furer 1997].
     * 
     * This version of the function is used as a subroutine in mod_exp_c_apx();
     */
    Set* domset = new Set();
    //TODO
    
    return domset;
}


//helpers
Set* max_card_sets(Graph* graph, Set* removedC, Set* removedS, int q) {
    /* Returns a set containing q vertices.
     * 
     * These vertices 
     * 
     */
    Set* q_verts = new Set();
    //TODO
    
    return q_verts;
}


bool p_cardinality(Graph* graph, Set* removedC, Set* removedS, int p) { 
    //Returns true IF all the residual subsets have cardinality at most p.
    bool run_apx=true;
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int v=*it;
        if(!removedS->contains(v) && !removedC->contains(v)) {
            Set* nbs_v = graph->neighbors(v);
            int sizeof_S=1; //1 for the vertex v itself
            
            for(auto inb=nbs_v->begin(); inb!=nbs_v->end(); inb++) {
                int nb_v=*inb;
                if(!removedC->contains(nb_v)) sizeof_S++;
            }
            
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
    
    for(auto it=graph->begin(); it!=graph->end(); it++) { //m*n
        int u=*it;
        
        if(!removedC->contains(u)) {
            int deg_u=0;
            Set* nbs_u = graph->neighbors(u);
            for(auto itt=nbs_u->begin(); itt!=nbs_u->end(); itt++) {
                int nb=*itt;
                if(!removedC->contains(nb)) deg_u++;     
            }
            
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
            int u_deg=0;
            for (auto in = graph->neighbors(u)->begin(); in!=graph->neighbors(u)->end(); in++) {
                int nb=*in;
                if(!visited->contains(nb)) u_deg++;
            }
            
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
