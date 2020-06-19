
#include "treewidth.hpp"

#include <fstream>
#include <limits>
#include <math.h>       // sqrt, log 
#include <cmath>        // floor


Set* treewidth_nodeedit(Graph* graph, int w) {
    /*
     * TreeWidthNodeEdit() from SR paper. 
     * 
     * rough psuedocode
     */
    
    int n = graph->size();
    
    int beta = floor(3*n/4);
    int c1 = 1;                 // constant value
    
    
    Set* Z = new Set();
    Set* W = new Set(); // not sure what to initialize these to.
    int t = tree_decomp(graph, Z, W);
    delete empty1, empty2;
    
    if ( t <= 32*c1*w*sqrt(log(w)) ) return 0 
        
    Set* S = balanced_separators(graph, beta);
    
    for componenets in G[V\S] {
         return S->set_union(treewidth_nodeedit(component, w)) ;
    }
}


int tree_decomp(Graph* graph, Set* Z, Set* W) {
    /*
     * Algorithm 4 from SR.
     * 
     * if Z ∩ W = ∅, output contains W in root bag.
     */
    int l = 0;
    
    if (8*Z->size() <= W->size()) {
        //return a tree decomposition with a single node containing Z ∩ W
    } else {
     
        int betaS = floor(3*W->size()/4);
        int betaT = floor(3*Z->set_union(W)->size()/4);
        Set* S = balanced_separators of W in G[Z ∪ W];
        Set* T = balanced_separators of Z ∪ W in G[Z ∪ W];
        
        l = num connected components;
        let G[V1], · · · , G[Vl] be the connected components of G[(W ∪ Z) \ (S ∪ T)
    }
    
    for (int i = 0; i < l; i++) {
        Set* Vi = connected component i;
        Set* Zi = Z->set_intersection(Vi);
        Set* Wi = W->set_intersection(Vi);
        
        Set* Ti = tree_decomp(graph, Zi, Wi->set_union(S)->set_union(T));
    }
    
    return tree decomposition with (W ∪ S ∪ T) as its root and T1, · · · , Tl as its children;
}


Set* balanced_separators(Graph* graph, int beta) {
    /*
     * Balanced separator greedy algorithm from from (Althoby et al. 2020)
     * 
     * Beta = floor(2n/3), or floor(3n/4)
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
        
        //Finds i in V\A  s.t.  |N(i) intersect B| is minimum.
        for (auto iu = graph->begin(); iu !=graph->end(); iu++) {
            int u = *iu;
            if(!A->contains(u)) {
                Set* nbi_B = graph->neighbors(u)->set_intersection(B);
                
                if (nbi_B->size() < min_inbsB)  {
                    min_inbsB = nbi_B->size();
                    vert = u;
                }
                delete nbi_B;
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
    
    printf("Sizeof A: %d\n", A->size());
    printf("Sizeof B: %d\n", B->size());
    printf("Sizeof C: %d\n", C->size());
    
    //1=true, 0=false
    printf("Separates? %d\n", test_separators(graph, A, B));
    
    delete A;
    delete B;
    
    return C; //Set of separator vertices
}


Set* balanced_separators(Graph graph, Set* W, int beta) {
    /*
     * Definition 5.1. For a subset of vertices W, a set of vertices S ⊆ V(G) is a vertex c-separator of W in G if each component of G[V \ S] contains at most c|W| vertices of W. The minimum size vertex c-separator of a graph, denoted sepc(G), is the minimum integer k such that for any subset W ⊆ V there exists a vertex c-separator of W in G of size k.
     * 
     */
    
    
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










