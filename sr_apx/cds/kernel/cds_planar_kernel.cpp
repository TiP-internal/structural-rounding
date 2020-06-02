
#include "cds_planar_kernel.hpp"
#include <vector>


Set* set_union(Set* set1, Set* set2) {
    //NOTE: Look into optimizing
     Set* set = new Set();
    
    for (Set::Iterator iu = set1->begin(); iu != set1->end(); ++iu) {
        int u = *iu;
        if (!set->contains(u)) set->insert(u);
    }
    
    for (Set::Iterator iu = set2->begin(); iu != set2->end(); ++iu) {
        int u = *iu;
        if (!set->contains(u)) set->insert(u);
    }
    
    return set;
}

Set* set_intersection(Set* set1, Set* set2) {
    //NOTE: Look into optimizing
    // set1 intersect set2
    
    Set* set = new Set();
    
    for (Set::Iterator iu = set1->begin(); iu != set1->end(); ++iu) {
        int u = *iu;
        if (set2->contains(u)) set->insert(u);
    }
    
    return set;
    
}

Set* set_minus(Set* set1, Set* set2) {
    //NOTE: Look into optimizing
    //  set1\set2
    Set* set = new Set();
    
    for (Set::Iterator iu = set1->begin(); iu != set1->end(); ++iu) {
        int u = *iu;
        if (!set2->contains(u)) set->insert(u);
    }
    
    return set;
}

Set* neighbors(Graph* graph, int u, int v) {
    //NOTE: Look into optimizing
    Set* pair_neighbors = new Set();
    
    for (Set::Iterator iv = graph->neighbors(v)->begin(); iv != graph->neighbors(v)->end(); ++iv) {
        int nv = *iv;
        if (nv != u && nv != v) {
            pair_neighbors->insert(nv);
        }
    }
    
    for (Set::Iterator iu = graph->neighbors(u)->begin(); iu != graph->neighbors(u)->end(); ++iu) {
        int nu = *iu;
        if (nu != u && nu != v) {
            pair_neighbors->insert(nu);
        }
    }
    
    return pair_neighbors;
}

Set* neighbors_closed(Graph* graph, int u) {
    //NOTE: Look into optimizing
    // need a deep copy here as to not modify the the adjacency list
    Set* closed = new Set();  
    for (Set::Iterator iu = graph->neighbors(u)->begin(); iu != graph->neighbors(u)->end(); ++iu) {
        closed->insert(*iu);
    }

    closed->insert(u);
    return closed;
}

Set* neighbors_closed(Graph* graph, int u, int v) {
    //NOTE: Look into optimizing
    Set* nbs = neighbors(graph, u, v);
    Set* closed = new Set();
    
    for (Set::Iterator nb = nbs->begin(); nb != nbs->end(); ++nb) {
        closed->insert(*nb);
    }
    
    closed->insert(u);
    closed->insert(v);
    return closed;
}


/*
 * Crust, mantle and core are the neighborhood partitioning functions used in [Gu & Imani, 2010].
 * They take the neighbohood (nbs) of a single vertex v and partition the neighbors 
 * into three partitions.
 */

Set* crust(Graph* graph, Set* nbs, Set* nbs_closed) {
    /*
     * nsb: the open neighborhood of vertex v
     * nbs_closed: closed neighborhood of v
     * 
     * Set of neighbors of v which have neighbors which are not also neighbors of v.
     */
    Set* crust_part = new Set();
    
    for (Set::Iterator ui = nbs->begin(); ui != nbs->end(); ++ui) {
        int u = *ui;
        
        Set* set = set_minus(graph->neighbors(u), nbs_closed);
        if (!set->empty()) crust_part->insert(u);
    }
    return crust_part;
}

Set* mantle(Graph* graph, Set* nbs, Set* crust_part) {
    /*
     * nbs: open neighborhood of vertex v
     * 
     * Set of neighbors of v which are also neighbors with vertices in the crust partition.
     */
    Set* mantle_part = new Set();
    Set* set_min = set_minus(nbs, crust_part);
    
    for (Set::Iterator ui = set_min->begin(); ui != set_min->end(); ++ui) {
        int u = *ui; 
        
        Set* set_inter = set_intersection(graph->neighbors(u), crust_part);
        if (!set_inter->empty()) mantle_part->insert(u);
    }
    return mantle_part;
}

Set* core(Set* nbs, Set* crust_part, Set* mantle_part) {
    /*
     * nbs: open neighborhood of vertex v
     * 
     * Set of neighbors of v which are not neighbors with the crust_part or the mantle_part.
     */
    Set* setun = set_union(crust_part, mantle_part);
    Set* core_part = set_minus(nbs, setun);
    return core_part;
}



