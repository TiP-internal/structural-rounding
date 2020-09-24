
#include "treewidth.hpp"
// #include "treedecomp.hpp"

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
    if(n==0) return components;
    
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


Set* treewidth_nodeedit(Graph* graph, Set* optional_verts, int w, bool annotated_version) {
    /*
     * TreeWidthNodeEdit() from SR paper. 
     * 
     */  
    int beta = floor(3*graph->size()/4);
    int c1 = 1;                         // NOTE constant value
    
    Set* W = new Set();
    Set* V = graph->get_vertices();
    
    TreeDecomp* decomp = new TreeDecomp();
    
    calculated_treedecomposition(graph, decomp);
    int t = decomp->treewidth();
    
    delete decomp;                    
    
    
    if ( t <= 32*c1*w*sqrt(log2(w)) ) {  //NOTE double check which log base
        Set* empty = new Set();
        return empty;
    }
    else { 
        Set* S = balanced_separators(graph, beta);  
        
        // for connected components of G[V\S]
        Set* V_min_S = V->set_minus(S);
        Graph* sub_g = graph->subgraph_wsingles(V_min_S);
        
        std::vector<Set*> components = connected_components(sub_g);
        delete V_min_S, sub_g;
        
        for (auto ic=components.begin(); ic!=components.end(); ic++) {
            Set* component_set = *ic;
            Graph* component = graph->subgraph_wsingles(component_set);
            
            Set* edited_vertices = treewidth_nodeedit(component, optional_verts, w, annotated_version);
            S->add_all(edited_vertices); //set union but modifies S
            
            if(annotated_version) {
                //Add the neighbors of the edit set to the optinally dominated set. 
                for(auto iev=edited_vertices->begin(); iev!=edited_vertices->end(); iev++) {
                    int e=*iev;
                    for(auto ien=component->neighbors(e)->begin(); 
                        ien!=component->neighbors(e)->end(); ien++) {
                        optional_verts->insert(*iev);
                    }
                }
            }
            
            delete component_set, component;
        }
        return S;
    }
}


void calculated_treedecomposition(Graph* graph, TreeDecomp* decomp) {
    // tree decomp doesnt quite work correctly when there are multiple components in graph.
    
    std::vector<Set*> components = connected_components(graph);
    int n_components = components.size();
    decomp->add_components(n_components);

    int i = 0;
    if(n_components > 1) {
        for (auto ic=components.begin(); ic!=components.end(); ic++) {
            
            Set* component_set = *ic;
            Graph* component = graph->subgraph_wsingles(component_set);
            
            Set* W = new Set();
            Set* V = component->get_vertices();
            
            tree_decomp(component, V, W, decomp->components_po_stacks[i], decomp->components_bags[i]);
            i++;
        }
    } else {
        Set* W = new Set();
        Set* V = graph->get_vertices();
        
        tree_decomp(graph, V, W, decomp->components_po_stacks[0], decomp->components_bags[0]);
    }
}


Set* tree_decomp(Graph* graph, Set* Z, Set* W, std::deque<std::deque<int>> &po_stack, std::vector<Set*> &bags) {
    /*
    * Algorithm 4 from SR.
    * 
    * Z is initialized to V.
    * W is initialized to empty set.
    * 
    * Bags stored in PRE-order traversal order.
    * 
    * ----------------------- NOTE: this is the bug fix! (though it's not very efficient yet)
    * 
    * The bug occurred when graph_Z_un_W formed a clique consisting of vertices which 
    * did not satisfy 8*Z->size() <= W->size(). 
    * 
    * The way the balanced separator were then calculated was by selecting one of the 
    * vertices as one of the 'shores,' while the rest were added the separator set, 
    * and the other 'shore' was empty. 
    * 
    * So |(S ∪ T)| == |graph_Z_un_W| - 1 and, 
    * this led to there being one connected component of G[(W ∪ Z) \ (S ∪ T)] (sub_g) consisting of one single 
    * vertex, which was the 'shore' found from the balanced separator (we'll call this vertex V). 
    * 
    * these sets were then input into the recursive call of tree_decomp, and the same 
    * separators and components were found, which led to an infinite number of recursive calls. 
    * 
    * NOTE this is fixed by essentially adding all vertices in the clique to the separator set, which
    * is done in the balanced separators function at the very end. 
    */ 
    
    std::vector<Set*> components;

    Set* T = new Set();
    Set* S_un_T = new Set();
    Set* Z_un_W = Z->set_union(W);
    
    int betaS, betaT;
       
    if(8*Z->size() <= W->size()) {  
        return Z_un_W;  
    } else {
        betaS = floor(3*W->size()/4);
        betaT = floor(3*Z_un_W->size()/4);
        
        Graph* graph_Z_un_W = graph->subgraph_wsingles(Z_un_W);
        
        //balanced_separators of W in G[Z ∪ W];
        //NOTE: balanced_separator calls take up something like 98% of the total time the alg. takes to run.
        Set* S = balanced_separators(graph_Z_un_W, W, betaS);  

        /* NOTE this is what takes the most time 
         * balanced_separators of Z ∪ W in G[Z ∪ W];
         * T = balanced_separators(graph_Z_un_W, Z_un_W, betaT);   
         * Switched to the version below. Less complicated function.
         */
        Set* T = balanced_separators(graph_Z_un_W, betaT);  
        
        // let G[V1], · · · , G[Vl] be the connected components of G[(W ∪ Z) \ (S ∪ T)]
        S_un_T = S->set_union(T);
        
        Set* Z_un_W_minus_S_un_T = Z_un_W->set_minus(S_un_T); 
        
        Graph* sub_g = graph->subgraph_wsingles(Z_un_W_minus_S_un_T);
        
        components = connected_components(sub_g);

        delete S, T, sub_g, Z_un_W_minus_S_un_T, graph_Z_un_W;
    }
    
    for (auto ic=components.begin(); ic!=components.end(); ic++) {
        Set* Vi = *ic;
    
        Set* Zi = Z->set_intersection(Vi);         
        Set* Wi = W->set_intersection(Vi);
        Set* Wi_un_S_un_T = Wi->set_union(S_un_T); 

        delete Vi, Wi;

        Set* Ti = tree_decomp(graph, Zi, Wi_un_S_un_T, po_stack, bags); 
        delete Zi, Wi_un_S_un_T;
        
        if(Ti->size()>0) {
            bags.push_back(Ti);
            int index = po_stack[0].size();               //Ti's parent's index
            
            if(index==0) { //no roots have been added yet
                //NOTE pre-order should NOT visit leaves first--figure out whats going on. 
                std::deque<int> leaves;
                po_stack.push_back(leaves); //empty deque, for the children of parent=bags.size()-1
            }
            
            po_stack[index].push_back(bags.size()-1);     //add Ti to the leaf of parent at index.
        } 
    }
    
    Set* W_un_S_un_T = W->set_union(S_un_T);  
    
    delete S_un_T, W, Z, Z_un_W;
    
    // return tree decomposition with (W ∪ S ∪ T) as its root and T1, · · · , Tl as its children;
    if(W_un_S_un_T->size()>0) { 
        bags.push_back(W_un_S_un_T);
        po_stack[0].push_back(bags.size()-1);  //add the parent to root array
        
        std::deque<int> leaves;
        po_stack.push_back(leaves); //empty deque, for the children of parent=bags.size()-1
        
        //add the childroot to leaves.
        int size = po_stack.size();

        if(size > 2) {  //add parent to end of leaves deque 
            po_stack[size-2].push_back(po_stack[0].back());
        }
    }
    
    Set* temp = new Set();
    return temp;
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
                nbi_B->remove(u);  //in case u is a single vertex in subgraph (neighbor of itself)

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
            if(*iv_nbs != vert) {    //in case vert is a single vertex in subgraph (neighbor of itself)
                if(!A->contains(*iv_nbs)) { 
                    C->insert(*iv_nbs); 
                    B->remove(*iv_nbs);
                }
            }
        }
    }
    
    //NOTE to double check the separators separate correctly.
    //bool test = test_separators(graph, A, B, beta);
    //printf("is valid sep? %d\n", test);
    
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
    int min_deg_v, min_deg_value;
    
    if(graph->size() == 1) { 
        // Single vertex must be the separator. Otherwise decomp doesnt work. 
        Set* separator = new Set();
        separator->insert(*graph->begin());
        return separator;
    }
    
    if(beta == 0 && W->size()==0) {  //no separator necessary
        Set* empty = new Set();
        return empty;
    }
    
    if(beta == 0 && W->size()>0) {  // beta is zero, then all of W must be in the set C (as part of separators).
        for (auto iw = W->begin(); iw!=W->end(); iw++) {
            int w = *iw;
            
            //find the min deg neighbor of w to add to A.
            float md_v = std::numeric_limits<float>::infinity();
            int min_deg_value = (int)md_v;
            int min_deg_vertex=0;
            
            //if W has no neighbors
            if(graph->neighbors(w)->size() <=1 && *graph->neighbors(w)->begin()==w) {
                if(!C->contains(w)) {
                    C->insert(w);
                    if(W->contains(w)) C_count++;
                }
            } else {  //else add min degree nebs of W to A, so W will be added to C.
                for (auto in=graph->neighbors(w)->begin(); in!=graph->neighbors(w)->end(); in++) {
                    int neb = *in;
                    if(neb != w) {
                        int neb_deg = graph->neighbors(neb)->size();
                        
                        if ( neb_deg < min_deg_value) {
                            min_deg_value = neb_deg;
                            min_deg_vertex = neb;
                        }
                    }
                }
                if(!A->contains(min_deg_vertex)) {
                    A->insert(min_deg_vertex);
                    if(W->contains(min_deg_vertex)) A_count++;
                }
            }
        }
    } else {  
        min_deg_v = min_deg_vert(graph);    
        A->insert(min_deg_v);
        if(W->contains(min_deg_v)) A_count++;
    }
    
    // C=N(A), where A={min_deg_v} or A=min deg neighbors of W.
    for (auto ia = A->begin(); ia != A->end(); ia++) {
        int a = *ia;
        for (auto ia_nbs=graph->neighbors(a)->begin(); ia_nbs!=graph->neighbors(a)->end(); ia_nbs++) {
            if(*ia_nbs != a) {  //in case min_deg_v is a single vertex in subgraph (neighbor of itself)
                if(!C->contains(*ia_nbs)) {
                    C->insert(*ia_nbs); 
                    if(W->contains(*ia_nbs)) C_count++;
                }
            }
        }
    }
    
     // B = V\(A union C)
    Set* AunC = A->set_union(C);
    for (auto ib = graph->begin(); ib!=graph->end(); ib++) {
        int b = *ib;
        if (!AunC->contains(b)) {
            if(!B->contains(b)) {
                B->insert(b);
                if (W->contains(b)) B_count++;
            }
        }
    }
    delete AunC;
    
    // if (W->size()==0), shouldnt enter while loop
    while(A_count + C_count < W->size() - beta) {        
        float min_inbsB = std::numeric_limits<float>::infinity();
        int vert; // the i s.t. |N(i) intersect B| is minimum.
        
        //Finds i in V\A  s.t.  |N(i) intersect B| is minimum.
        for (auto iu = graph->begin(); iu !=graph->end(); iu++) {
            int u = *iu;
            if(!A->contains(u)) {
                Set* nbi_B = graph->neighbors(u)->set_intersection(B);
                nbi_B->remove(u);  //in case u is a single vertex in subgraph (neighbor of itself)
                
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
        
        for (auto iv_nbs = graph->neighbors(vert)->begin();
                iv_nbs != graph->neighbors(vert)->end(); iv_nbs++) {
            
            if(*iv_nbs != vert) {  //in case vert is a single vertex in subgraph (neighbor of itself)
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
    
    
    if(B->size() == 0) {
        //NOTE figure out the size of graph for which to test this. As to minimize the number of times this is called
        bool is_graph_clique = is_clique(graph);   //NOTE: not the most efficient function either.
        if(is_graph_clique) {    //NOTE: in case of future issues, potentially get rid of this case. 
            for (auto ia = A->begin(); ia!=A->end(); ia++) {
                int a = *ia;
                C->insert(a);
                
                A->remove(a);
                if (W->contains(a)) A_count--;
            }
        }
    }
    
    //NOTE For testing. Tests whether the vertex separators are valid separators.
    //bool test = test_separators(graph, A, B, A_count, B_count, beta);
    //if (!test) {

    //    printf("W: ");
    //    for(auto it=W->begin(); it!=W->end(); it++) printf(" %d,", *it);
    //    printf("\n\n");
        
    //    printf("beta=%d\n", beta);
    //    printf("A_count=%d, B_count=%d, C_count=%d \n", A_count, B_count, C_count);
    //    printf("ERROR: is valid sep? %d\n", test);
        
    //}
    
    delete A;
    delete B;
    
    return C; //Set of separator vertices
}



/*      Testing Functions --now necessary for balanced seps      */

bool test_separators(Graph* graph, Set* A, Set* B, int A_count, int B_count, int beta) {
    /*
     * Tests that there are not edges between the two sets.
     */
    //printf("is |A|=%d <=? beta=%d\n", A->size(), beta);
    //printf("is |B|=%d <=? beta=%d\n", B->size(), beta);
    
    //if(A->size() > beta) return false;  
    //if(B->size() > beta) return false;
    
    if(A_count > beta) {
        printf("A: ");
        for(auto it=A->begin(); it!=A->end(); it++) printf(" %d,", *it);
        printf("\n\n");
        
        return false;
    }
    if(B_count > beta) {
        return false;
    }
    
    for (Set::Iterator iu = A->begin(); iu != A->end(); ++iu) {
        int u = *iu;

        for (Set::Iterator iv = B->begin(); iv != B->end(); ++iv) {
            int v = *iv;
            if (graph->adjacent(u,v)) {
                printf("EDGE\n");
                return false; //there is an edge between the sets
            }
        }
    }
    return true;
}


bool is_clique(Graph* graph) {
 
    bool clique = true;
    
    for (auto iu = graph->begin(); iu !=graph->end(); iu++) {
        int u = *iu;
        
        for (auto iv = graph->begin(); iv !=graph->end(); iv++) {
            int v = *iv;
            
            if(u!= v && !graph->adjacent(u,v)) clique = false;
        }
    }
    
    return clique;
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






