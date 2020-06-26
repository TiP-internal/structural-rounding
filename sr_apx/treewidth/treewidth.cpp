
#include "treewidth.hpp"

#include <fstream>
#include <limits>
#include <math.h>       // sqrt, log 
#include <cmath>        // floor
#include <vector>


int* treewidth_nodeedit(Graph* graph, Set* edit_set, Set* Z, Set* W, std::vector<Set*> &bags, int w) {
    /*
     * TreeWidthNodeEdit() from SR paper. 
     * 
     */

    int n = graph->size();    
    int beta = floor(3*n/4);
    
    printf("n: %d, beta:%d\n", n, beta);
    for(auto iv=graph->begin(); iv!=graph->end(); iv++) {
            printf("%d, ", *iv);
    }
    printf("\n");
    
    int c1 = 1;                 // constant value
    
    //bags.push_back(W);
    tree_decomp(graph, Z, W, bags);
    int t = 3;
    
    
    //if ( t <= 32*c1*w*sqrt(log(w)) ) return post_order();
    //printf("c1*w*sqrt(log(w))=%f\n", c1*w*sqrt(log(w)));
    if ( n < 1 ) return post_order();
        
    Set* S = balanced_separators(graph, beta);
    printf("Size of S: %d\n", S->size());
    for (auto ie=S->begin(); ie!=S->end(); ie++) {
        printf("     S v: %d\n", *ie);
    }
    
    edit_set->add_all(S) ;    //set union but modifies S
    
    for (auto ie=edit_set->begin(); ie!=edit_set->end(); ie++) {
        printf("     edit_set v: %d\n", *ie);
    }
    
    // for connected components of G[V\S]
    std::vector<Set*> components = connected_components(graph, S);
    printf("number of components: %ld\n", components.size());
    for (auto ic=components.begin(); ic!=components.end(); ic++) {
        Set* c = *ic;
        for (auto iv=c->begin(); iv!=c->end(); iv++) {
            printf("  v: %d, ", *iv);
        }
        printf("\n_____\n");
    }
    
    for (auto ic=components.begin(); ic!=components.end(); ic++) {
        Set* component_set = *ic;
        Graph* component = graph->subgraph(component_set);
        
        treewidth_nodeedit(component, edit_set, W, Z, bags, w);
        delete component;
    }
    delete S;
    
    return 0;
}


std::vector<Set*> connected_components(Graph* graph, Set* S) {
    /*
     * DFS for finding connected components in G[V\S]
     */
    
    int n = graph->size();
    std::vector<Set*> components;
    Set* visited = new Set();
    auto vitr = graph->begin();
    
    while(vitr!=graph->end()) {  
        if(!visited->contains(*vitr) && !S->contains(*vitr)) {  //if not visited yet
            Set* comp = new Set();
            dfs(graph, comp, visited, S, *vitr);
            components.push_back(comp);
        }
        ++vitr;
    }
    delete visited;
    return components;
}

void dfs(Graph* graph, Set* comp, Set* visited, Set* S, int v) {
    visited->insert(v);

    if (!S->contains(v)) comp->insert(v);  // components in G[V\S]
    
    Set* nbrs = graph->neighbors(v);
    for (auto it = nbrs->begin(); it != nbrs->end(); it++) {
        if(!visited->contains(*it) && !S->contains(*it)) {    //if not visited yet
            dfs(graph, comp, visited, S, *it);
        }
    }
}


Set* tree_decomp(Graph* graph, Set* Z, Set* W, std::vector<Set*> &bags) {
    /*
     * Algorithm 4 from SR.
     * 
     * if Z ∩ W = ∅, output contains W in root bag.
     * 
     * Z is initialized to empty set.
     * W is initialized to empty set, but should be the V.
     * 
     * But maybe actually its the opposite? Z=V, and W=empty
     */
    printf("\n---------------------------tree_Decomp\n");
    
    int Z_size;
    
    //NOTE: make sure that this is the only time Z is empty.
    if(Z->size()==0) {  // first iteration, Z=V (but Z is initialized to empty to avoid an additional loop)
        Z_size = graph->size();
        printf("W_size=%d, Z_size=%d\n", W->size(), Z_size);
        if(8*Z_size <= W->size()) { //NOTE double check that this happens?--if not get rid of if stmt
            printf("First return\n");
            printf("---------------------------\n\n");
            //return a tree decomposition with a single node containing Z ∩ W
            Set* empty_set = new Set(); //Z should be empty here.
            //bags.push_back(empty_set);
            //return -1; 
            return empty_set;
        }
    } else{
        Z_size = Z->size();
        printf("W_size=%d, Z_size=%d\n", W->size(), Z_size);
        if(8*Z_size <= W->size()) { 
            //return a tree decomposition with a single node containing Z ∩ W
            Set* Z_int_W = Z->set_intersection(W);
            
            printf("Second return\n");
            printf("Sizeof Z_int_W=%d\n", Z_int_W->size());
            bags.push_back(Z_int_W);
            printf("---------------------------\n\n");
            
            //int tw = Z_int_W->size()-1;            
            return Z_int_W; 
        }
    }
    
    Set* Z_un_W = Z->set_union(W);
    int betaS = floor(3*W->size()/4);
    int betaT = floor(3*Z_un_W->size()/4);
    
    printf("Z_un_W size=%d\n", Z_un_W->size());
    printf("betaS=%d, betaT=%d\n",betaS, betaT);
        
    Graph* graph_Z_un_W = graph->subgraph(Z_un_W);
    
    //balanced_separators of W in G[Z ∪ W];
    Set* S = balanced_separators(graph_Z_un_W, W, betaS);
    
    //balanced_separators of Z ∪ W in G[Z ∪ W];
    Set* T = balanced_separators(graph_Z_un_W, Z_un_W, betaT);

    delete Z_un_W;
    
    // let G[V1], · · · , G[Vl] be the connected components of G[(W ∪ Z) \ (S ∪ T)]
    Set* S_un_T = S->set_union(T);
    Set* Z_un_W_minus_S_un_T = Z_un_W->set_minus(S_un_T);
    
    std::vector<Set*> components = connected_components(graph, Z_un_W_minus_S_un_T);
    for (auto ic=components.begin(); ic!=components.end(); ic++) {
        Set* Vi = *ic;

        Set* Zi = Z->set_intersection(Vi);
        Set* Wi = W->set_intersection(Vi);
        
        Set* Wi_un_S_un_T = Wi->set_union(S_un_T);        
        Set* Ti = tree_decomp(graph, Zi, Wi_un_S_un_T, bags);
        
        printf("Sizeof Ti=%d\n", Ti->size());
        bags.push_back(Ti);
    }
    delete S_un_T, Z_un_W_minus_S_un_T;
    // return tree decomposition with (W ∪ S ∪ T) as its root and T1, · · · , Tl as its children;

    Set* temp = new Set();
    printf("******Returning temp\n");
    return temp;
}


int* post_order() {
    printf("returning post_order()\n");
    return 0;
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
    
//     printf("Sizeof A: %d\n", A->size());
//     printf("Sizeof B: %d\n", B->size());
//     printf("Sizeof C: %d\n", C->size());
//     
//     //1=true, 0=false
//     printf("Separates? %d\n", test_separators(graph, A, B));
    
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
    
    printf("**bal seps: %d, beta: %d\n", W->size(), beta);
    
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
        C->insert(*ia_nbs); 
        if(W->contains(*ia_nbs)) C_count++;
    }
         
         
    // B = V\(A union C)
    Set* AunC = A->set_union(C);
    for (auto ib = graph->begin(); ib != graph->end(); ib++) {
        int b = *ib;
        if (!AunC->contains(b)) {
            B->insert(b);
            if (W->contains(b)) B_count++;
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
        A->insert(vert);
        if (W->contains(vert)) A_count++;
        
        B->remove(vert);
        if (W->contains(vert)) B_count--;
        
        // C=N(A), where A = A union vert
        C->remove(vert);
        if(W->contains(vert)) C_count--;
        for (Set::Iterator iv_nbs = graph->neighbors(vert)->begin();
            iv_nbs != graph->neighbors(vert)->end(); iv_nbs++) {
            
            if(!A->contains(*iv_nbs)) { 
                C->insert(*iv_nbs); 
                B->remove(*iv_nbs);
                if (W->contains(*iv_nbs)) B_count--;
            }
        }
    }
    
//     printf("A_Count=%d, B_Count=%d, C_Count=%d\n", A_count, B_count, C_count);
//     
//     printf("\n\nSizeof A: %d\n", A->size());
//     
//     for(auto iv=A->begin(); iv!=A->end(); iv++) {
//         printf("   A_v: %d\n", *iv);
//     }
//     
//     printf("\n\nSizeof B: %d\n", B->size());
//     
//     for(auto iv=B->begin(); iv!=B->end(); iv++) {
//         printf("   B_v: %d\n", *iv);
//     }
//     
//     printf("\n\nSizeof C: %d\n", C->size());
//     
//     for(auto iv=C->begin(); iv!=C->end(); iv++) {
//         printf("   C_v: %d\n", *iv);
//     }
//     
//     //1=true, 0=false
//     printf("Separates? %d\n", test_separators(graph, A, B));
    
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










