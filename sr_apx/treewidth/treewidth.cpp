
#include "treewidth.hpp"
#include "treedecomp.hpp"

#include "time.h"

#include <fstream>
#include <limits>
#include <math.h>       // sqrt, log 
#include <cmath>        // floor
#include <vector>


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
    nbrs->remove(v);
    for (auto it = nbrs->begin(); it != nbrs->end(); it++) {
        if(!visited->contains(*it)) {    //if not visited yet
            dfs(graph, comp, visited, *it);
        }
    }
}


Set* treewidth_nodeedit(Graph* graph, int w) {
    /*
     * TreeWidthNodeEdit() from SR paper. 
     * 
     */  
    int beta = floor(3*graph->size()/4);
    int c1 = 1;                         // NOTE constant value
    
    Set* W = new Set();
    Set* V = graph->get_vertices();
    
    TreeDecomp* decomp = new TreeDecomp();
    
    tree_decomp(graph, V, W, decomp->preorder_stack, decomp->bags);
    int t = decomp->treewidth();
    delete decomp;                      // NOTE delete or store?
    
    
    if ( t <= 32*c1*w*sqrt(log2(w)) ) {  //NOTE double check which log base
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

bool sets_equal(Set* A, Set* B) {
    bool are_equal = true;
    
    for(auto iv=A->begin(); iv!=A->end(); iv++) {
        if(!B->contains(*iv)) are_equal=false;
    }

    for(auto iv=B->begin(); iv!=B->end(); iv++) {
        if(!A->contains(*iv)) are_equal=false;
    }
    
    return are_equal;
}


Set* tree_decomp(Graph* graph, Set* Z, Set* W, std::deque<std::deque<int>> &preorder_stack, std::vector<Set*> &bags) {
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
    
    if(graph->size() == 1) {  //in case a component of size one is input.
        Set* single_bag = new Set();
        single_bag->insert(*graph->begin());
        bags.push_back(single_bag);
        preorder_stack[0].push_back(bags.size()-1);
    }
    
    if(8*Z->size() <= W->size()) {  
        return Z_un_W;  
    } else {
        int betaS = floor(3*W->size()/4);
        int betaT = floor(3*Z_un_W->size()/4);
        
        Graph* graph_Z_un_W = graph->subgraph(Z_un_W);
        
        //NOTE: balanced_separator calls take up something like 98% of the total time the alg. takes to run.
        //balanced_separators of W in G[Z ∪ W];
        S = balanced_separators(graph_Z_un_W, W, betaS);  

        //balanced_separators of Z ∪ W in G[Z ∪ W];
        T = balanced_separators(graph_Z_un_W, Z_un_W, betaT);   //NOTE this is what takes the most time
        
        // let G[V1], · · · , G[Vl] be the connected components of G[(W ∪ Z) \ (S ∪ T)]
        S_un_T = S->set_union(T);
        
        Set* Z_un_W_minus_S_un_T = Z_un_W->set_minus(S_un_T); 
        
        Graph* sub_g = graph->subgraph(Z_un_W_minus_S_un_T);
        components = connected_components(sub_g);
        
        delete S, T, sub_g, graph_Z_un_W, Z_un_W_minus_S_un_T;
    }
    
    int Z_size = Z->size();
    int W_size = W->size();
    int n_components = components.size();
    std::deque<int> leaves;
    for (auto ic=components.begin(); ic!=components.end(); ic++) {
        Set* Vi = *ic;
    
        Set* Zi = Z->set_intersection(Vi);         
        Set* Wi = W->set_intersection(Vi);
        delete Vi;  // deleting this causes probelsm
        
        Set* Wi_un_S_un_T = Wi->set_union(S_un_T); 
        
        //-----NOTE potential fix for the infinite loop problem. Doesnt guarantee that |Z| <= |W|/8 though.
        if(sets_equal(Z, Zi) && sets_equal(W, Wi_un_S_un_T)) {
            return Z_un_W;
        }
        //-----
        
        Set* Ti = tree_decomp(graph, Zi, Wi_un_S_un_T, preorder_stack, bags); 
        
        delete Zi, Wi, Wi_un_S_un_T;
        
        if(Ti->size()>0) {
            bags.push_back(Ti);
            leaves.push_back(bags.size()-1);
        }
    }
    
    preorder_stack.push_back(leaves);
    Set* W_un_S_un_T = W->set_union(S_un_T);  
    
    // return tree decomposition with (W ∪ S ∪ T) as its root and T1, · · · , Tl as its children;
    if(W_un_S_un_T->size()>0) {
        bags.push_back(W_un_S_un_T);
        preorder_stack[0].push_back(bags.size()-1);
    }
    
    delete S_un_T, W, Z, Z_un_W;
    //NOTE: need to return something here, but is this it?
    Set* temp = new Set();
    return temp;
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
    float md_v = std::numeric_limits<float>::infinity();
    int min_deg_value = (int)md_v;
    int min_deg_vert=0;
    
    //Find a minimum degree vertex to add to A
    for (auto iu = graph->begin(); iu!=graph->end(); iu++) {
        int u = *iu;
        int u_deg = graph->neighbors(u)->size();
        
        if ( u_deg < min_deg_value) {
            min_deg_value = u_deg;
            min_deg_vert = u;
        }
    }
    return min_deg_vert;
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
    
    int min_deg_v = min_deg_vert(graph);;

    A->insert(min_deg_v);
    
    // C=N(A), where A={min_deg_v}
    for (Set::Iterator ia_nbs = graph->neighbors(min_deg_v)->begin();
         ia_nbs != graph->neighbors(min_deg_v)->end(); ia_nbs++) {
        if(*ia_nbs != min_deg_v) C->insert(*ia_nbs); 
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
                nbi_B->remove(u);

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
            if(*iv_nbs != vert) {
                if(!A->contains(*iv_nbs)) { 
                    C->insert(*iv_nbs); 
                    B->remove(*iv_nbs);
                }
            }
        }
    }
//     bool test = test_separators(graph, A, B, beta);
//     printf("is valid sep? %d\n", test);
    
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
    int min_deg_v;
    
    if (W->size() == 1) {
        //here W must be in the separator, floor(3/4)=0, meaning each component can have 0 elements of W. 
        //shouldnt enter while loop below
        min_deg_v = *W->begin();
    } else {
        min_deg_v = min_deg_vert(graph);
    }
    
    A->insert(min_deg_v);
    if(W->contains(min_deg_v)) A_count++;
    
//     for(auto iv=A->begin(); iv!=A->end(); iv++) printf("A: %d, ", *iv);
//     printf("\n");
    
    // C=N(A), where A={min_deg_v}
    for (Set::Iterator ia_nbs = graph->neighbors(min_deg_v)->begin(); ia_nbs != graph->neighbors(min_deg_v)->end(); ia_nbs++) {
        if(*ia_nbs != min_deg_v) {
            if(!C->contains(*ia_nbs)) {
                C->insert(*ia_nbs); 
                if(W->contains(*ia_nbs)) C_count++;
            }
        }
    }
    
//     for(auto iv=C->begin(); iv!=C->end(); iv++) printf("C: %d, ", *iv);
//     printf("\n");
            
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
    
//     for(auto iv=B->begin(); iv!=B->end(); iv++) printf("B: %d, ", *iv);
//     printf("\n");
    
//     printf("**A_count=%d, C_count=%d, B_count=%d\n", A_count, C_count, B_count);
//     printf("**W->size=%d, beta=%d\n", W->size(), beta);
//     printf("**A_count + C_count=%d, <? W->size() - beta= %d\n", (A_count + C_count), (W->size() - beta));
    
    
    // if (W->size()==0), shouldnt enter while loop
    while(A_count + C_count < W->size() - beta) { 
//         printf("    A_count=%d, C_count=%d, B_count=%d\n", A_count, C_count, B_count);
//         printf("    W->size=%d, beta=%d\n", W->size(), beta);
//         printf("    A_count + C_count=%d, <? W->size() - beta= %d\n", (A_count + C_count), (W->size() - beta));
        
        float min_inbsB = std::numeric_limits<float>::infinity();
        int vert; // the i s.t. |N(i) intersect B| is minimum.
        
        //Finds i in V\A  s.t.  |N(i) intersect B| is minimum.
        for (auto iu = graph->begin(); iu !=graph->end(); iu++) {
            int u = *iu;
            if(!A->contains(u)) {
                Set* nbi_B = graph->neighbors(u)->set_intersection(B);
                nbi_B->remove(u);
                
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
            
            if(*iv_nbs != vert) {
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
    }
    
//     bool test = test_separators(graph, A, B, beta);
//     printf("is valid sep? %d\n", test);
    
    delete A;
    delete B;
    
    return C; //Set of separator vertices
}


/*      Testing Functions       */

bool test_separators(Graph* graph, Set* A, Set* B, int beta) {
    /*
     * Tests that there are not edges between the two sets.
     */
    printf("is |A|=%d <=? beta=%d\n", A->size(), beta);
    printf("is |B|=%d <=? beta=%d\n", B->size(), beta);
    
    for (Set::Iterator iu = A->begin(); iu != A->end(); ++iu) {
        int u = *iu;

        for (Set::Iterator iv = B->begin(); iv != B->end(); ++iv) {
            int v = *iv;
            if (graph->adjacent(u,v)) return false; //there is an edge between the sets
        }
    }
    return true;
}





