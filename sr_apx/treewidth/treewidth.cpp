
#include "treewidth.hpp"
// #include "treedecomp.hpp"

#include "time.h"

#include <fstream>
#include <limits>
#include <math.h>       // sqrt, log
#include <cmath>        // floor
#include <vector>

TreeDecomp::TreeDecomp() {
    tw = -1;
}

TreeDecomp::TreeDecomp(Graph* graph): TreeDecomp() {
    build_decomposition(graph);
}

TreeDecomp::~TreeDecomp() {
    for (Set* s : components_bags) {
        delete s;
    }
}

int TreeDecomp::treewidth() {
    return tw - 1;
}

std::vector<po_bag> TreeDecomp::get_post_order() {
    std::vector<po_bag> po;
    for (auto it = pre_order.rbegin(); it != pre_order.rend(); ++it) {
        po.push_back(*it);
    }

    return po;
}

Set* copyset(Set* s) {
    Set* copy = new Set();
    for (int x : *s) {
        copy->insert(x);
    }
    return copy;
}

// returns placement of new bag
// parent_index should be 0 for root bag
// only modifies pre_order and components_bags
int TreeDecomp::add_bag(int parent, bool last_child, Set* bag) {
    // ensures that root bag exists
    if (components_bags.empty()) {
        components_bags.push_back(new Set());
        pre_order.push_back(po_bag(0, -1));
    }

    int true_parent = pre_order[parent].current_join_child;

    if (pre_order[true_parent].num_children > 0 && !last_child) {
        // create a new join bag to hold the child
        int index = components_bags.size();
        components_bags.push_back(copyset(components_bags[true_parent]));
        pre_order.push_back(po_bag(index, true_parent));
        pre_order[true_parent].num_children++;
        pre_order[parent].current_join_child = index;
        true_parent = index;
    }

    if (pre_order[true_parent].num_children > 0 || !last_child) {
        int index = components_bags.size();
        components_bags.push_back(copyset(components_bags[true_parent]));
        pre_order.push_back(po_bag(index, true_parent));
        pre_order[true_parent].num_children++;
        true_parent = index;
    }

    // create the child
    int use_parent = true_parent;
    Set* use_bag = copyset(components_bags[true_parent]);

    // add introduce bags
    for (int x : *(components_bags[true_parent])) {
        // only look at vertices to be added
        if (bag->contains(x)) {
            continue;
        }

        use_bag->erase(x);
        int index = components_bags.size();
        components_bags.push_back(copyset(use_bag));
        pre_order.push_back(po_bag(index, use_parent));
        pre_order[use_parent].num_children++;
        use_parent = index;
    }

    // add forget bags
    for (int x : *bag) {
        // only look at vertices to be removed
        if (components_bags[true_parent]->contains(x)) {
            continue;
        }

        use_bag->insert(x);
        int index = components_bags.size();
        components_bags.push_back(copyset(use_bag));
        pre_order.push_back(po_bag(index, use_parent));
        pre_order[use_parent].num_children++;
        use_parent = index;
    }

    tw = tw < bag->size() ? bag->size() : tw;
    delete use_bag;

    // return index of final bag
    return components_bags.size() - 1;
}

void TreeDecomp::build_decomposition(Graph* graph) {
    // tree decomp doesnt quite work correctly when there are multiple components in graph.

    std::vector<Set*> components = connected_components(graph);
    int n_components = components.size();


    if(n_components > 1) {
        for (int i = 0; i + 1 < n_components; ++i) {
            Set* component_set = components[i];
            Graph* component = graph->subgraph_wsingles(component_set);

            Set* W = new Set();
            Set* V = component->get_vertices();

            tree_decomp(component, V, W, 0, false);
            delete component_set;
            delete component;
            delete V;
            delete W;
        }

        Set* component_set = components[n_components - 1];
        Graph* component = graph->subgraph_wsingles(component_set);

        Set* W = new Set();
        Set* V = component->get_vertices();

        tree_decomp(component, V, W, 0, true);
        delete component_set;
        delete component;
        delete W;
        delete V;

    }
    else {
        Set* W = new Set();
        Set* V = graph->get_vertices();

        tree_decomp(graph, V, W, 0, true);
        delete components[0];
        delete V;
        delete W;
    }
}

void TreeDecomp::tree_decomp(Graph* graph, Set* Z, Set* W, int parent, bool last_child) {
    /*
    * Algorithm 4 from SR.
    *
    * Z is initialized to V.
    * W is initialized to empty set.
    *
    * Bags stored in PRE-order traversal order.
    *
    */

    std::vector<Set*> components;

    Set* Z_un_W = Z->set_union(W);

    int betaS, betaT;

    if(8*Z->size() <= W->size()) {
        int leaf = add_bag(parent, last_child, Z_un_W);
        delete Z_un_W;
        return;
    }

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
    Set* S_un_T = S->set_union(T);

    Set* Z_un_W_minus_S_un_T = Z_un_W->set_minus(S_un_T);

    Graph* sub_g = graph->subgraph_wsingles(Z_un_W_minus_S_un_T);

    components = connected_components(sub_g);

    delete S;
    delete T;
    delete sub_g;
    delete Z_un_W_minus_S_un_T;
    delete graph_Z_un_W;

    Set* W_un_S_un_T = W->set_union(S_un_T);
    int x = add_bag(parent, last_child, W_un_S_un_T);
    delete W_un_S_un_T;

    for (int i = 0; i + 1 < components.size(); i++) {
        Set* Vi = components[i];

        Set* Zi = Z->set_intersection(Vi);
        Set* Wi = W->set_intersection(Vi);
        Set* Wi_un_S_un_T = Wi->set_union(S_un_T);

        delete Vi;
        delete Wi;
        tree_decomp(graph, Zi, Wi_un_S_un_T, x, false);
        delete Zi;
        delete Wi_un_S_un_T;
    }

    if (components.size() > 0) {
        int i = components.size() - 1;
        Set* Vi = components[i];

        Set* Zi = Z->set_intersection(Vi);
        Set* Wi = W->set_intersection(Vi);
        Set* Wi_un_S_un_T = Wi->set_union(S_un_T);

        delete Vi;
        delete Wi;
        tree_decomp(graph, Zi, Wi_un_S_un_T, x, true);
        delete Zi;
        delete Wi_un_S_un_T;
    }

    delete S_un_T;
    delete Z_un_W;
}

// end of treedecomp object //////////////////////////

std::vector<Set*> connected_components(Graph* graph) {
    /*
     * DFS for finding connected components in G[V]
     *
     * n+m time
     */

    std::vector<Set*> components;

    std::vector<int> stack;
    Set visited;
    Map<Set>::Iterator vitr = graph->begin();

    Set* component = NULL;

    while (!stack.empty() || visited.size() < graph->size()) {
        int current;
        if (stack.empty()) {
            if (component != NULL) {
                components.push_back(component);
            }

            component = new Set();

            while (visited.contains(*vitr)) {
                ++vitr;
            }

            current = *vitr;
            component->insert(current);
            visited.insert(current);
        }
        else {
            current = stack.back();
            stack.pop_back();
        }

        for (int nbr : *(graph->neighbors(current))) {
            if (!component->contains(nbr)) {
                stack.push_back(nbr);
                component->insert(nbr);
                visited.insert(nbr);
            }
        }
    }

    if (component != NULL) {
        components.push_back(component);
    }

    return components;

}


Set* treewidth_nodeedit(Graph* graph, Set* optional_verts, int w, bool annotated_version) {
    /*
     * TreeWidthNodeEdit() from SR paper.
     *
     */
    int beta = floor(3*graph->size()/4);
    int c1 = 1;                         // NOTE constant value

    Set* V = graph->get_vertices();

    TreeDecomp decomp(graph);
    int t = decomp.treewidth();


    if ( t <= c1*w*sqrt(log2(w)) ) {  //NOTE double check which log base
        Set* empty = new Set();
        delete V;
        return empty;
    }

    Set* S = balanced_separators(graph, beta);

    // for connected components of G[V\S]
    Set* V_min_S = V->set_minus(S);
    delete V;
    Graph* sub_g = graph->subgraph_wsingles(V_min_S);

    std::vector<Set*> components = connected_components(sub_g);
    delete V_min_S;
    delete sub_g;

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
                    int optional = *ien;
                    if(!S->contains(optional)) {
                        optional_verts->insert(*ien);
                    }
                    optional_verts->remove(e);
                }
            }
        }

        delete edited_vertices;
        delete component_set;
        delete component;
    }

    return S;
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
        delete A;
        delete B;
        delete C;
        return separator;
    }

    if(beta == 0 && W->size()==0) {  //no separator necessary
        Set* empty = new Set();
        delete A;
        delete B;
        delete C;
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
        //NOTE figure out the size of graph for which to test this. 
        //As to minimize the number of times this is called
        //NOTE: not the most efficient function either.
        bool is_graph_clique = is_clique(graph);   
        if(is_graph_clique) {    
            //NOTE: in case of future issues, potentially get rid of this case.
            for (auto ia = A->begin(); ia!=A->end(); ia++) {
                int a = *ia;
                C->insert(a);

                A->remove(a);
                if (W->contains(a)) A_count--;
            }
        }
    }

    delete A;
    delete B;

    return C; //Set of separator vertices
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
