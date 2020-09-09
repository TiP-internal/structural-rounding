
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
     */
    Set* domset;
    Set* visited = new Set();
    
    int p=largest_harmonic_num(q);
    while(visited->size() < graph->size()) {
        
        //1. IF there exists an item of C that belongs to a single subset S ∈ S, THEN add S to the solution
        
        
        //2. IF there exist two sets S, R in S such that S is included into R, THEN remove S without branching;
        
        
        /* 3. IF all the residual subsets have cardinality at most p, 
         *    THEN run the algorithm by [Duh & Furer '97] in order to compute a q-approximation of the 
         *    optimal solution in the surviving instance
         */  
        domset = Hk_minus_half_apx(graph, visited);
        
        /* 4. determine q sets S1, . . . , Sq from S such that (union i<=q Si) has maximum cardinality and perform 
         *    the the following branching: 
         * 
         *      a. either add every Si to the solution (and remove ∪i?qSi from C), 
         *      b. or remove all of them.
         */ 
        
    }
    
    delete visited;
    return domset;
}


Set* Hk_minus_half_apx(Graph* graph, Set* visited) {  
    /*
     * Approximation algorithm (originally for set cover) from "Approximationof k-SetCover by Semi-Local Optimization"
     * [Duh & Furer 1997].
     * 
     * This version of the function is used as a subroutine in mod_exp_c_apx();
     */
    Set* domset = new Set();
    
    
    return domset;
}


//helpers
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
