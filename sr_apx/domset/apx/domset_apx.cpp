
#include <cstdio>

#include "domset_apx.hpp"

int max_deg_vertex(Graph* graph, Set* visited) {
    //finds the maximum degree not in the visited set
    int max_deg_value = 0;
    int max_deg_vert=0;
    
    for (auto iu = graph->begin(); iu!=graph->end(); iu++) {
        int u = *iu;
        if(!visited->contains(u)) {
            int u_deg = graph->neighbors(u)->size();
            
            if ( u_deg > max_deg_value) {
                max_deg_value = u_deg;
                max_deg_vert = u;
            }
        }
    }
    return max_deg_vert;
}


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
    
    return domset;
}
