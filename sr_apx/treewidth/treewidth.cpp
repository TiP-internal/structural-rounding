
#include "treewidth.hpp"
#include "treedecomp.hpp"

#include "time.h"

#include <fstream>
#include <limits>
#include <math.h>       // sqrt, log 
#include <cmath>        // floor
#include <vector>


Set* treewidth_nodeedit(Graph* graph, int w) {
    /*
     * TreeWidthNodeEdit() from SR paper. 
     * 
     */
    int n = graph->size();    
    int beta = floor(3*n/4);
    int c1 = 1;                 // NOTE constant value
    
    Set* W = new Set();
    Set* V = graph->get_vertices();
    
    TreeDecomp* decomp = new TreeDecomp();
    
    tree_decomp(graph, V, W, decomp->preorder_stack, decomp->bags, 0);
    int t = decomp->treewidth();
    delete decomp;              // NOTE delete or store?
    
//     if (t <= w) {   // NOTE for testing
//         Set* empty = new Set();
//         return empty;
//     } 
    if ( t <= 32*c1*w*sqrt(log(w)) ) {
        Set* empty = new Set();
        return empty;
    }
    else { 
        Set* S = balanced_separators(graph, beta);    
        
        // for connected components of G[V\S]
        Set* V_min_S = V->set_minus(S);
        Graph* sub_g = graph->subgraph(V_min_S);
        
        std::vector<Set*> components = connected_components(sub_g);
        delete V_min_S, sub_g;
        
        for (auto ic=components.begin(); ic!=components.end(); ic++) {
            Set* component_set = *ic;
            Graph* component = graph->subgraph(component_set);
            
            S->add_all(treewidth_nodeedit(component, w)); //set union but modifies S

            delete component_set, component;
        }
        return S;
    }
}


std::vector<Set*> connected_components(Graph* graph) {
    /*
     * DFS for finding connected components in G[V]
     * 
     * n+m time
     */
    
    int n = graph->size();
    std::vector<Set*> components;
    Set* visited = new Set();
    auto vitr = graph->begin();
    
    while(vitr!=graph->end()) {  
        if(!visited->contains(*vitr)) {  //if not visited yet
            Set* comp = new Set();
            dfs(graph, comp, visited, *vitr);
            components.push_back(comp);
        }
        ++vitr;
    }
    delete visited;
    return components;
}

void dfs(Graph* graph, Set* comp, Set* visited, int v) {
    visited->insert(v);
    comp->insert(v);  // components in G[V]
    
    Set* nbrs = graph->neighbors(v);
    for (auto it = nbrs->begin(); it != nbrs->end(); it++) {
        if(!visited->contains(*it)) {    //if not visited yet
            dfs(graph, comp, visited, *it);
        }
    }
}

Set* tree_decomp(Graph* graph, Set* Z, Set* W, std::deque<std::deque<int>> &preorder_stack, std::vector<Set*> &bags, int po_index) {
    /*
     * Algorithm 4 from SR.
     * 
     * Z is initialized to V.
     * W is initialized to empty set.
     * 
     * Bags stored in POST-order traversal order.
     * 
     * NOTE this doesnt work when graph has more than one connected components
     */
    
    std::vector<Set*> components;
    
    Set* S = new Set();
    Set* T = new Set();
    Set* S_un_T = new Set();
    Set* Z_un_W = Z->set_union(W);
    
    //WARNING this is the problem im talking about. 8*Z->size is greater than W->size, but Z doesnt get
    // reduced in size by the code before the recursive call, resulting in infinte recursive calls. 
    printf("8*Z->size()=%d, W->size()=%d\n", 8*Z->size(), W->size());
    
    if(8*Z->size() <= W->size()) {
//     if(Z->size() <= W->size()) {        //NOTE for testing only
        return Z_un_W;
    } else {
        int betaS = floor(3*W->size()/4);
        int betaT = floor(3*Z_un_W->size()/4);
            
        Graph* graph_Z_un_W = graph->subgraph(Z_un_W);
        
        //balanced_separators of W in G[Z ∪ W];
        S = balanced_separators(graph_Z_un_W, W, betaS);
        
        //balanced_separators of Z ∪ W in G[Z ∪ W];
        T = balanced_separators(graph_Z_un_W, Z_un_W, betaT);
        
        // let G[V1], · · · , G[Vl] be the connected components of G[(W ∪ Z) \ (S ∪ T)]
        S_un_T = S->set_union(T);
        Set* Z_un_W_minus_S_un_T = Z_un_W->set_minus(S_un_T);  //NOTE could replace w/ remove_all fun.
        
        Graph* sub_g = graph->subgraph(Z_un_W_minus_S_un_T);
        components = connected_components(sub_g);
        
        delete S, T, Z_un_W, Z_un_W_minus_S_un_T, sub_g, graph_Z_un_W;
    }
    
    int index = 0;
    int n_components = components.size();
    std::deque<int> leaves;  
    for (auto ic=components.begin(); ic!=components.end(); ic++) {
        Set* Vi = *ic;
    
        Set* Zi = Z->set_intersection(Vi);   //WARNING Here, Zi is not reduced in size. //NOTE could also replace w/ same_elements fun.
        
        Set* Wi = W->set_intersection(Vi);
        delete Vi;
        
        Set* Wi_un_S_un_T = Wi->set_union(S_un_T);  //NOTE replace w. add_all fun.      
        
        Set* Ti = tree_decomp(graph, Zi, Wi_un_S_un_T, preorder_stack, bags, po_index++);  //NOTE recursive call
        delete Zi;
        
        if(Ti->size()>0) {
            bags.push_back(Ti);
            leaves.push_back(po_index-1);
        }
        delete Wi, Wi_un_S_un_T;
        index++;
    }
    
    preorder_stack.push_back(leaves);
    Set* W_un_S_un_T = W->set_union(S_un_T);   //NOTE replace w. add_all fun.
    delete S_un_T;

    // return tree decomposition with (W ∪ S ∪ T) as its root and T1, · · · , Tl as its children;
    preorder_stack[0].push_back(po_index-n_components);
    bags.push_back(W_un_S_un_T);

    
    //NOTE: need to return something here, but is this it?
    Set* temp = new Set();
    return temp;
    //return W_un_S_un_T;  // Doesnt work
}

int find_treewidth(std::vector<Set*> &bags) {
    int tw = 0;
    for(auto ib=bags.begin(); ib!=bags.end(); ib++) {
            Set* bag = *ib;
            if (bag->size() > tw) tw = bag->size();
    }
    return tw-1;
}

int min_deg_vert(Graph* graph) {
    float min_deg_v = std::numeric_limits<float>::infinity();
    
    //Find a minimum degree vertex to add to A
    for (auto iu = graph->begin(); iu!=graph->end(); iu++) {
        int u = *iu;
        if ( u < min_deg_v) {
            min_deg_v = u;
        }
    }
    return min_deg_v;
}


Set* balanced_separators(Graph* graph, int beta) {
    /*
     * NOTE: could probably get rid of this function, and replace call 
     * on line 31 w/ other balanced_separators(graph, V, beta);
     * 
     * Balanced separator greedy algorithm from from (Althoby et al. 2020)
     * 
     * Beta = floor(2n/3), or floor(3n/4)
     * 
     */
        
    Set* A = new Set();
    Set* B = new Set();
    Set* C = new Set();
    
    float min_deg_v = min_deg_vert(graph); 
    
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
    delete A;
    delete B;
    
    return C; //Set of separator vertices
}


Set* balanced_separators(Graph* graph, Set* W, int beta) {
    /*
     * Balanced separator greedy algorithm from from (Althoby et al. 2020)
     * modified to find separators of subset W. If W=V it's the same. 
     * 
     * Beta = floor(2n/3), or floor(3n/4)
     * 
     * Definition 5.1. For a subset of vertices W, a set of vertices S ⊆ V(G) 
     * is a vertex c-separator of W in G if each component of G[V \ S] contains 
     * at most c|W| vertices of W. The minimum size vertex c-separator of a graph, 
     * denoted sepc(G), is the minimum integer k such that for any subset W ⊆ V 
     * there exists a vertex c-separator of W in G of size k.
     * 
     */
    int A_count=0, B_count=0, C_count=0;  // counts of vertices in A,B,C which are in W.
    Set* A = new Set();
    Set* B = new Set();
    Set* C = new Set();
    
    float min_deg_v = min_deg_vert(graph);
    
    A->insert(min_deg_v);
    if(W->contains(min_deg_v)) A_count++;
    
    // C=N(A), where A={min_deg_v}
    for (Set::Iterator ia_nbs = graph->neighbors(min_deg_v)->begin();
         ia_nbs != graph->neighbors(min_deg_v)->end(); ia_nbs++) {
        if(!C->contains(*ia_nbs)) {
            C->insert(*ia_nbs); 
            if(W->contains(*ia_nbs)) C_count++;
        }
    }
         
         
    // B = V\(A union C)
    Set* AunC = A->set_union(C);
    for (auto ib = graph->begin(); ib != graph->end(); ib++) {
        int b = *ib;
        if (!AunC->contains(b)) {
            if(!B->contains(b)) {
                B->insert(b);
                if (W->contains(b)) B_count++;
            }
        }
    }
    delete AunC;
    
    while(A_count + C_count < W->size() - beta) { //N
        
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
        if(!A->contains(vert)) {
            A->insert(vert);
            if (W->contains(vert)) A_count++;
        }
        
        if (B->contains(vert)) {
            B->remove(vert);
            if (W->contains(vert)) B_count--;
        }
        
        // C=N(A), where A = A union vert
        if(C->contains(vert)) {
            C->remove(vert);
            if(W->contains(vert)) C_count--;
        }
        
        for (Set::Iterator iv_nbs = graph->neighbors(vert)->begin();
            iv_nbs != graph->neighbors(vert)->end(); iv_nbs++) {
            
            if(!A->contains(*iv_nbs)) { 
                if (!C->contains(*iv_nbs)) {
                    C->insert(*iv_nbs);
                    if(W->contains(*iv_nbs)) C_count++;
                }
                
                if (B->contains(*iv_nbs)) {
                    B->remove(*iv_nbs);
                    if (W->contains(*iv_nbs)) B_count--;
                }
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










