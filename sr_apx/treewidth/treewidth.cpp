
#include "treewidth.hpp"

#include <fstream>
#include <limits>


Set* balanced_separators(Graph* graph, int beta) {
    /*
     * Balanced separator greedy algorithm from from (Althoby et al. 2020)
     * 
     * Beta = floor(2n/3)
     */
        
    Set* A = new Set();
    Set* B = new Set();
    Set* C = new Set();
    
    float min_deg_v = std::numeric_limits<float>::infinity();
    
    //Find a minimum degree vertex to add to A
    for (auto iu = graph->begin(); iu!=graph->end(); iu++) {
        int u = *iu;
        if ( u < min_deg_v) {
            min_deg_v = u;
        }
    }
    
    A->insert(min_deg_v);
    
    // C=N(A), where A={min_deg_v}
    for (Set::Iterator ia_nbs = graph->neighbors(min_deg_v)->begin();
         ia_nbs != graph->neighbors(min_deg_v)->end(); ia_nbs++) {
        C->insert(*ia_nbs); 
    }
         
         
    // B = V\(A union C)
    Set* AunC = A->set_union(C);
    for (auto ib = graph->begin(); ib != graph->end(); ib++) {
        int b = *ib;
        if (!AunC->contains(b)) B->insert(b);
    }
    delete AunC;
    
    
    while(A->size() + C->size() < graph->size() - beta) {
        float min_inbsB = std::numeric_limits<float>::infinity();
        int vert; // the i s.t. |N(i) intersect B| is minimum.
        
        //Finds Let i in V\A  s.t.  |N(i) intersect B| is minimum.
        for (auto iu = graph->begin(); iu !=graph->end(); iu++) {
            int u = *iu;
            if(!A->contains(u)) {
                Set* nbi_B = graph->neighbors(u)->set_intersection(B);
                
                if (nbi_B->size() < min_inbsB)  {
                    min_inbsB = nbi_B->size();
                    vert = u;
                }
                delete nbi_B;  //necessary?
            }
        }
        
        //A = A union vert
        A->insert(vert);
        B->remove(vert);
        
        // C=N(A), where A = A union vert
        C->remove(vert);
        for (Set::Iterator iv_nbs = graph->neighbors(vert)->begin();
            iv_nbs != graph->neighbors(vert)->end(); iv_nbs++) {
            
            if(!A->contains(*iv_nbs)) { 
                C->insert(*iv_nbs); 
                B->remove(*iv_nbs);
            }
        }
        
    }    
    delete A;
    delete B;
    
    return C; //Set of separator vertices
}


bool test_separators(Graph* graph, Set* A, Set* B) {
    /*
     * Tests that there are not edges between the two sets.
     */
    for (Set::Iterator iu = A->begin(); iu != A->end(); ++iu) {
        int u = *iu;
        
        for (Set::Iterator iv = B->begin(); iv != B->end(); ++iv) {
            int v = *iv;
            if (graph->adjacent(u,v)) return false; //there is an edge between the sets
        }
    }
    return true;
}










