
#include "sr_apx/treewidth/treewidth.hpp"

namespace sr_apx {
namespace treewidth {

Decomposition::Decomposition(bool b = true) {
    width = 0;
    build = b;
}

Decomposition::Decomposition(const Graph& graph): Decomposition(true) {
    build_decomposition(graph);
}

int Decomposition::treewidth() {
    return width - 1;
}

// returns placement of new bag
// parent_index should be 0 for root bag
// only modifies pre_order and components_bags
int Decomposition::add_bag(int parent, bool last_child, const Set& bag) {
    width = width < bag.size() ? bag.size() : width;

    if (!build) {
        return -1;
    }

    // ensures that root bag exists
    if(pre_order.empty()) {
        pre_order.push_back(po_bag());
    }

    int true_parent = pre_order[parent].current_join_child;

    if (pre_order[true_parent].left_child != -1 && !last_child) {
        // create a new join bag to hold the child
        int index = pre_order.size();

        pre_order.push_back(po_bag(index, pre_order[true_parent].bag));
        pre_order[parent].current_join_child = index;
        pre_order[true_parent].right_child = index;
        true_parent = index;
    }

    // create the child
    int index = pre_order.size();
    std::vector<int> temp_bag;
    for (int b : bag) {
        temp_bag.push_back(b);
    }
    pre_order.push_back(po_bag(index, std::move(temp_bag)));

    if (pre_order[true_parent].left_child == -1) {
        pre_order[true_parent].left_child = index;
    }
    else {
        pre_order[true_parent].right_child = index;
    }

    // return index of final bag
    return index;
}

void Decomposition::build_decomposition(const Graph& graph) {
    build_decomposition(graph, graph.size() + 1);
}

void Decomposition::build_decomposition(const Graph& graph, int limit) {
    // std::vector<Set> components = graph.connected_components();
    // int nc = components.size();
    //
    // for (int i = 0; i + 1 < nc; ++i) {
    //     Graph component = graph.subgraph(components[i]);
    //     Set empty;
    //     tree_decomp(component, components[i], empty, 0, false);
    // }
    //
    // Graph component = graph.subgraph(components[nc - 1]);
    // Set empty;
    // tree_decomp(component, components[nc - 1], empty, 0, true);

    Map<int> first_bag;
    std::vector<int> stack;

    Map<Set>::const_iterator vertices = graph.begin();

    while (first_bag.size() < graph.size()) {
        // early out for decision variant
        if (width > limit + 1) {
            return;
        }

        Graph subgraph;
        Set W;
        int parent_bag = 0;

        for ( ; vertices != graph.end(); ++vertices) {
            if (!first_bag.contains(vertices->first)) {
                stack.push_back(vertices->first);
                break;
            }
        }

        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();

            if (first_bag.contains(u)) {
                W.insert(u);
                parent_bag = parent_bag > first_bag[u] ? parent_bag : first_bag[u];
                continue;
            }

            for (int v : graph.neighbors(u)) {
                if (!subgraph.contains_vertex(v)) {
                    stack.push_back(v);
                }
                subgraph.add_edge(u, v);
            }
        }

        if (subgraph.size() == 0) {
            Set S;
            S.insert(vertices->first);
            int index = add_bag(parent_bag, false, S);
            first_bag[vertices->first] = index;
            continue;
        }

        if (subgraph.size() <= width) {
            Set S;
            Map<Set>::const_iterator iu = subgraph.begin();
            for ( ; iu != subgraph.end(); ++iu) {
                S.insert(iu->first);
            }

            int index = add_bag(parent_bag, false, S);
            for (int s : S) {
                if (!first_bag.contains(s)) {
                    first_bag[s] = index;
                }
            }
            continue;
        }

        Set sep = balanced_separator(subgraph, W);
        sep.insert(W.begin(), W.end());
        int index = add_bag(parent_bag, false, sep);
        for (int s : sep) {
            if (!first_bag.contains(s)) {
                first_bag[s] = index;
            }
        }
    }
}

void Decomposition::build_gd_decomposition(const Graph& graph) {
    build_gd_decomposition(graph, graph.size() + 1);
}

void Decomposition::build_gd_decomposition(const Graph& graph, int limit) {
    /* Algorithm 3 from Treewidth computations I. Upper bounds */

    Graph h = graph; /* h = graph */
    int n = h.size();

    /* initialize deg_sets */
    std::vector<Set> deg_sets(n, Set());
    for (Map<Set>::const_iterator it = h.begin(); it != h.end(); ++it) {
        int degree = it->second.size();
        deg_sets[degree].insert(it->first);
    }

    /* iterate over h */
    int early_out = 0;
    int min_degree = 1;
    std::vector<int> ordering;
    std::vector<Set> bags;
    ordering.reserve(n);
    bags.reserve(n);
    for (int i = 0; i < n; i++) {
        if (early_out > limit + 1) {
            width = early_out;
            return; /* early out for decision variant */
        }

        int v = min_vertex(h, deg_sets, min_degree); /* choose v as a vertex of smallest degree in h */
        min_degree = h.degree(v);
        early_out = min_degree + 1 > early_out ? min_degree + 1 : early_out;
        if (min_degree < 1) min_degree = 1;
        if (!h.contains_vertex(v)) { /* hotfix for graphs with self-loop nodes */
            i--;
            continue;
        }

        /* make a clique of v's neighbors in h */
        std::vector<int> neighbors;
        for (int w : h.neighbors(v))
            neighbors.push_back(w);

        for (int j = 0; j+1 < neighbors.size(); j++) {
            int w = neighbors[j];
            int w_deg = h.degree(w);
            int count = 0;
            for (int k = j+1; k < neighbors.size(); k++) {
                int x = neighbors[k];
                int x_deg = h.degree(x);
                if (!h.adjacent(w,x)) {
                    h.add_edge(w,x);
                    /* increase degree of x */
                    deg_sets[x_deg].erase(x);
                    deg_sets[x_deg+1].insert(x);
                    count++;
                }
            }
            /* increase degree of w */
            if (count > 0) {
                deg_sets[w_deg].erase(w);
                deg_sets[w_deg+count].insert(w);
            }
        }

        Set b = h.neighbors(v);
        b.insert(v);
        bags.push_back(b);

        /* remove v and all its edges */
        for (int w : h.neighbors(v)) {
            int degree = h.degree(w);
            deg_sets[degree].erase(w);
            deg_sets[degree-1].insert(w);
        }
        h.remove_vertex(v);

        ordering.push_back(v); /* let v be the last vertex in the ordering */
    }

    Map<int> bag_index;

    int index = add_bag(0, false, bags[n-1]);
    bag_index[ordering[n-1]] = index;

    for (int i = n - 2; i >= 0; --i) {
        int parent = 0;
        for (int nbr : bags[i]) {
            if (nbr == ordering[i]) {
                continue;
            }

            int x = bag_index[nbr];
            parent = x > parent ? x : parent;
        }

        index = add_bag(parent, false, bags[i]);
        bag_index[ordering[i]] = index;
    }
}

void Decomposition::build_gfi_decomposition(const Graph& graph) {
    build_gfi_decomposition(graph, graph.size()+1);
}

void Decomposition::build_gfi_decomposition(const Graph& graph, int limit) {
    /* Algorithm 3 from Treewidth Computations I. Upper bounds */
  
    Graph h = graph;
    int n = h.size();

    /* calculate max possible fill for any v */
    int max_fill = 1;
    for (int i = n-1; i > 0; i--)
        max_fill += i;

    /* initialize fill_sets */
    Map<int> fill(h.size());
    std::vector<Set> fill_sets(max_fill, Set());
    for (Map<Set>::const_iterator it = h.begin(); it != h.end(); ++it) {
        int f = fill_edges(h, it->first);
        fill_sets[f].insert(it->first);
        fill[it->first] = f;
    }

    /* iterator over h */
    int early_out = 0;
    int min_fill = 1;
    std::vector<int> ordering;
    std::vector<Set> bags;
    ordering.reserve(n);
    bags.reserve(n);
    for (int i = 0; i < n; i++) {
      //if (early_out > limit + 1) {
          //width = early_out;
          //return; /* early out for decision variant */
      //}

        int v = min_vertex(h, fill_sets, min_fill); /* choose v as a vertex that causes the smallest number of fill edges */
        //min_fill =;
        //early_out =;
        if (min_fill < 0) min_fill = 0;
        if (!h.contains_vertex(v)) { /* hotfix for graphs with self-loop nodes */
            i--;
            continue;
        }

        /* make a clique of v's neighbors in h */
        std::vector<int> neighbors;
        for (int w : h.neighbors(v))
            neighbors.push_back(w);

        for (int j = 0; j+1 < neighbors.size(); j++) {
            int w = neighbors[j];
            int count = 0;
            for (int k = j+1; k < neighbors.size(); k++) {
	        int x = neighbors[k];
	        if (!h.adjacent(w,x)) {
	            h.add_edge(w,x);
	            /* decrease the fill of shared neighbors */  
		    Set w_neighbors = h.neighbors(w);
		    Set x_neighbors = h.neighbors(x);
		    for (int y : w_neighbors) {
		        if (x_neighbors.contains(y) && y != v) {
			    fill_sets[fill[y]].erase(y);
			    int new_fill = fill[y]-1;
			    fill_sets[new_fill].insert(y);
			    fill[y] = new_fill;
			}
		    }
		}
            }
	}

        Set b = h.neighbors(v);
        b.insert(v);
        bags.push_back(b);

        /* remove v and all its edges */
        h.remove_vertex(v);
	for (int w : neighbors) {
	    fill_sets[fill[w]].erase(w);
	    int new_fill = fill_edges(h, w);
	    //printf("\nold fill[w] = %d\n", fill[w]);
	    //printf("new fill[w] = %d\n", new_fill);
	    fill_sets[new_fill].insert(w);
	    fill[w] = new_fill;
        }
        h.remove_vertex(v);

        ordering.push_back(v); /* let v be the last vertex in the ordering */
    }

    Map<int> bag_index;

    int index = add_bag(0, false, bags[n-1]);
    bag_index[ordering[n-1]] = index;

    for (int i = n - 2; i >= 0; --i) {
        int parent = 0;
        for (int nbr : bags[i]) {
            if (nbr == ordering[i]) {
	        continue;
            }

            int x = bag_index[nbr];
            parent = x > parent ? x : parent;
        }
    
        index = add_bag(parent, false, bags[i]);
        bag_index[ordering[i]] = index;
    }
}

void Decomposition::build_decomposition(const Graph& graph, std::vector<int> ordering) {
    Graph h = graph;
    tree_decomp_ordering(h, ordering.size(), ordering);
}

void Decomposition::tree_decomp_ordering(Graph& graph, int n, std::vector<int> ordering) {
    Graph g = graph;
    Set bags [n];

    /* add vn's bag containing only vn */
    Set x_vn;
    int vn = ordering[n-1];
    x_vn.insert(vn);
    g.remove_vertex(vn);
    bags[n-1] = x_vn;
    add_bag(0, false, x_vn);

    /* add the rest of the bags */
    for (int i = n-2; i > -1; i--) {
        Set x_vi;
        int vi = ordering[i];
        x_vi.insert(vi); /* add vi to bag */

        /* add vi's forward neighbors */
        for (int w : g.neighbors(vi))
            x_vi.insert(w);
        g.remove_vertex(vi);

        bags[i] = x_vi;
        if (i == 0) add_bag(ordering[i+1], true, x_vi);
        else add_bag(ordering[i+1], false, x_vi);
    }
}

void Decomposition::tree_decomp(Graph& graph, Set& Z, Set& W, int parent, bool last_child) {
    /*
    * Algorithm 4 from SR.
    *
    * Z is initialized to V.
    * W is initialized to empty set.
    *
    * Bags stored in PRE-order traversal order.
    *
    */

    Set Z_union_W = Z;
    Z_union_W.insert(W.begin(), W.end());

    if (Z_union_W.size() <= width) {
        add_bag(parent, last_child, Z_union_W);
        return;
    }

    Graph graph_Z_union_W = graph.subgraph(Z_union_W);
    Set sep = balanced_separator(graph_Z_union_W, W);

    W.insert(sep.begin(), sep.end());
    for (int s : sep) {
        Z.erase(s);
    }

    for (int s : sep) {
        for (int w : W) {
            graph.remove_edge(s, w);
        }
    }

    int x = add_bag(parent, last_child, W);

    if (Z.empty()) {
        return;
    }

    Graph sub_g = graph.subgraph(Z);
    std::vector<Set> components = sub_g.connected_components();
    int nc = components.size();

    for (int i = 0; i + 1 < nc; ++i) {
        Set Wi;
        for (int v : components[i]) {
            for (int nbr : graph.neighbors(v)) {
                if (W.contains(nbr)) {
                    Wi.insert(nbr);
                }
            }
        }
        tree_decomp(graph, components[i], Wi, x, false);
    }

    Set Wi;
    for (int v : components[nc - 1]) {
        for (int nbr : graph.neighbors(v)) {
            if (W.contains(nbr)) {
                Wi.insert(nbr);
            }
        }
    }
    tree_decomp(graph, components[nc - 1], Wi, x, true);
}

void Decomposition::sort_bags() {
    for (po_bag& bag : pre_order) {
        std::sort(bag.bag.begin(), bag.bag.end());
    }
}

// end of Decomposition object //////////////////////////

Set vertex_delete(const Graph& graph, int w) {
    /*
     * TreeWidthNodeEdit() from SR paper.
     *
     */

    std::vector<Set> comps = graph.connected_components();
    int nc = comps.size();

    if (nc > 1) {
        Set res;
        for (int i = 0; i < nc; i++) {
            Graph component = graph.subgraph(comps[i]);
            Set part = vertex_delete(component, w);
            res.insert(part.begin(), part.end());
        }

        return res;
    }

    Decomposition decomp(false);
    decomp.build_decomposition(graph, w);
    // decomp.build_gd_decomposition(graph, w);
    int t = decomp.treewidth();

    if (t <= w) {
        return Set();
    }

    Set S = balanced_separator(graph, Set());

    Set V_minus_S;
    Map<Set>::const_iterator iu = graph.begin();
    for ( ; iu != graph.end(); ++iu) {
        if (!S.contains(iu->first)) {
            V_minus_S.insert(iu->first);
        }
    }

    Graph sub_g = graph.subgraph(V_minus_S);
    Set res = vertex_delete(sub_g, w);
    res.insert(S.begin(), S.end());
    return res;
}

Set vertex_gd_delete(const Graph& graph, int w) {
    /*
     * TreeWidthNodeEdit() from SR paper.
     *
     */

    std::vector<Set> comps = graph.connected_components();
    int nc = comps.size();

    if (nc > 1) {
        Set res;
        for (int i = 0; i < nc; i++) {
            Graph component = graph.subgraph(comps[i]);
            Set part = vertex_gd_delete(component, w);
            res.insert(part.begin(), part.end());
        }

        return res;
    }

    Decomposition decomp(false);
    // decomp.build_decomposition(graph, w);
    decomp.build_gd_decomposition(graph, w);
    int t = decomp.treewidth();

    if (t <= w) {
        return Set();
    }

    Set S = balanced_separator(graph, Set());

    Set V_minus_S;
    Map<Set>::const_iterator iu = graph.begin();
    for ( ; iu != graph.end(); ++iu) {
        if (!S.contains(iu->first)) {
            V_minus_S.insert(iu->first);
        }
    }

    Graph sub_g = graph.subgraph(V_minus_S);
    Set res = vertex_gd_delete(sub_g, w);
    res.insert(S.begin(), S.end());
    return res;
}

Set vertex_gfi_delete(const Graph& graph, int w) {
    std::vector<Set> comps = graph.connected_components();
    int nc = comps.size();

    if (nc > 1) {
        Set res;
        for (int i = 0; i < nc; i++) {
            Graph component = graph.subgraph(comps[i]);
            Set part = vertex_gfi_delete(component, w);
            res.insert(part.begin(), part.end());
        }

        return res;
    }

    Decomposition decomp(false);
    decomp.build_gfi_decomposition(graph, w);
    int t = decomp.treewidth();

    if (t <= w) {
        return Set();
    }

    Set S = balanced_separator(graph, Set());

    Set V_minus_S;
    Map<Set>::const_iterator iu = graph.begin();
    for ( ; iu != graph.end(); ++iu) {
        if (!S.contains(iu->first)) {
            V_minus_S.insert(iu->first);
        }
    }

    Graph sub_g = graph.subgraph(V_minus_S);
    Set res = vertex_gfi_delete(sub_g, w);
    res.insert(S.begin(), S.end());
    return res;
}

Set balanced_separator(const Graph& graph, const Set& W) {
    if (graph.size() <= 2) {
        Set C;
        Map<Set>::const_iterator iu = graph.begin();
        for ( ; iu != graph.end(); ++iu) {
            C.insert(iu->first);
        }
        return C;
    }

    // cannot legitimately call with empty W so used for V case
    if (W.empty()) {
        Set V;
        Map<Set>::const_iterator iu = graph.begin();
        for ( ; iu != graph.end(); ++iu) {
            V.insert(iu->first);
        }
        return balanced_separator(graph, V);
    }

    if (W.size() < 2) {
        Set C;
        int w = *(W.begin());
        int nbr = *(graph.neighbors(w).begin());
        C.insert(nbr);
        return C;
    }

    Map<int> deg;
    std::vector<Set> revdeg;
    int maxdeg = 0;

    Map<Set>::const_iterator iu = graph.begin();
    for ( ; iu != graph.end(); ++iu) {
        // deg is number of neighbors in B
        int degree = iu->second.size() + 1;
        deg[iu->first] = degree;
        maxdeg = degree > maxdeg ? degree : maxdeg;
        revdeg.resize(maxdeg+1);
        revdeg[degree].insert(iu->first);
    }

    int A_count = 0;
    int C_count = 0;
    Set C;

    while (A_count + C_count < W.size() - ((2 * W.size()) / 3) || (!revdeg[0].empty() && A_count + C_count < W.size() - 1) || (!revdeg[1].empty() && A_count + C_count < W.size() - 1)) {
        int mindeg;
        for (mindeg = 0; revdeg[mindeg].empty(); ++mindeg);
        int vert = *(revdeg[mindeg].begin());

        // move vert to A
        revdeg[mindeg].erase(vert);
        deg.erase(vert);
        if (W.contains(vert)) {
            ++A_count;
            if (C.contains(vert)) {
                --C_count;
            }
        }

        for (int add_c : graph.neighbors(vert)) {
            // neighbor already in A
            if (!deg.contains(add_c)) {
                continue;
            }

            // update degree of add_c if vert was in B
            if (!C.contains(vert)) {
                int degree = deg[add_c];
                deg[add_c] = degree - 1;
                revdeg[degree].erase(add_c);
                revdeg[degree - 1].insert(add_c);
            }

            // don't readd to C
            if (C.contains(add_c)) {
                continue;
            }

            for (int nbr : graph.neighbors(add_c)) {
                // nbr in A
                if (!deg.contains(nbr)) {
                    continue;
                }

                int degree = deg[nbr];
                deg[nbr] = degree - 1;
                revdeg[degree].erase(nbr);
                revdeg[degree - 1].insert(nbr);
            }

            if (W.contains(add_c)) {
                ++C_count;
            }
            C.insert(add_c);

            int degree = deg[add_c];
            deg[add_c] = degree - 1;
            revdeg[degree].erase(add_c);
            revdeg[degree - 1].insert(add_c);
        }

        C.erase(vert);
    }

    return C;
}


Set close_balanced_separator(const Graph& graph, const Set& W) {
    if (graph.size() <= 2) {
        Set C;
        Map<Set>::const_iterator iu = graph.begin();
        for ( ; iu != graph.end(); ++iu) {
            C.insert(iu->first);
        }
        return C;
    }

    // cannot legitimately call with empty W so used for V case
    if (W.empty()) {
        Set V;
        Map<Set>::const_iterator iu = graph.begin();
        for ( ; iu != graph.end(); ++iu) {
            V.insert(iu->first);
        }
        return close_balanced_separator(graph, V);
    }

    if (W.size() < 2) {
        Set C;
        int w = *(W.begin());
        int nbr = *(graph.neighbors(w).begin());
        C.insert(nbr);
        return C;
    }

    Map<int> deg;
    std::vector<Set> revdeg;
    int maxdeg = 0;
    int mindeg = graph.size() + 1;
    int min_deg_vertex = -1;

    Map<Set>::const_iterator iu = graph.begin();
    for ( ; iu != graph.end(); ++iu) {
        int degree = iu->second.size() + 1;
        maxdeg = degree > maxdeg ? degree : maxdeg;
        if (degree < mindeg) {
            mindeg = degree;
            min_deg_vertex = iu->first;
        }
    }

    deg[min_deg_vertex] = 0;
    revdeg.resize(maxdeg + 1);
    revdeg[0].insert(min_deg_vertex);

    for (int nbr : graph.neighbors(min_deg_vertex)) {
        int degree = graph.degree(nbr) + 1;
        deg[nbr] = degree;
        revdeg[degree].insert(nbr);
    }

    int A_count = 0;
    int C_count = 0;
    Set C;
    Set A;

    while (A_count + C_count < W.size() - ((2 * W.size()) / 3) || (!revdeg[0].empty() && A_count + C_count < W.size() - 1) || (!revdeg[1].empty() && A_count + C_count < W.size() - 1)) {
        int mindeg;
        for (mindeg = 0; revdeg[mindeg].empty(); ++mindeg);
        int vert = *(revdeg[mindeg].begin());

        // move vert to A
        revdeg[mindeg].erase(vert);
        deg.erase(vert);
        A.insert(vert);
        if (W.contains(vert)) {
            ++A_count;
            if (C.contains(vert)) {
                --C_count;
            }
        }

        for (int add_c : graph.neighbors(vert)) {
            // neighbor already in A
            if (A.contains(add_c)) {
                continue;
            }

            if (!deg.contains(add_c)) {
                int degree = graph.degree(add_c) + 1;
                deg[add_c] = degree;
                revdeg[degree].insert(add_c);
            }

            // update degree of add_c if vert was in B
            if (!C.contains(vert)) {
                int degree = deg[add_c];
                deg[add_c] = degree - 1;
                revdeg[degree].erase(add_c);
                revdeg[degree - 1].insert(add_c);
            }

            // don't readd to C
            if (C.contains(add_c)) {
                continue;
            }

            for (int nbr : graph.neighbors(add_c)) {
                // nbr in A
                if (A.contains(nbr)) {
                    continue;
                }

                // nbr is a new available choice
                if (!deg.contains(nbr)) {
                    int degree = graph.degree(nbr) + 1;
                    deg[nbr] = degree;
                    revdeg[degree].insert(nbr);
                }

                int degree = deg[nbr];
                deg[nbr] = degree - 1;
                revdeg[degree].erase(nbr);
                revdeg[degree - 1].insert(nbr);
            }

            if (W.contains(add_c)) {
                ++C_count;
            }
            C.insert(add_c);

            int degree = deg[add_c];
            deg[add_c] = degree - 1;
            revdeg[degree].erase(add_c);
            revdeg[degree - 1].insert(add_c);
        }

        C.erase(vert);
    }

    return C;
}

// greedy heuristics

int chordal(Graph&, std::vector<Set>&);
// int chordal(const Graph&, std::vector<Set>&);
bool clique(const Graph&, int);
int lowest_neighbor(const Graph&, int, std::vector<int>);
int min_vertex(const Graph&, std::vector<Set>&, int);
int min_degree_vertex1(const Graph&);
int min_fill_edges_vertex1(const Graph&);


Graph fill(const Graph& g, int n, std::vector<int> ordering) {
    /* Algorithm 1 from Treewidth computations I. Upper bounds */

    /* create map of ordering */
    Map<int> pi;
    for (int i = 0; i < n; i++)
        pi[ordering[i]] = i;

    Graph h = g; /* h = g */

    for (int i = 0; i < n; i++) {
        int v = ordering[i]; /* let v be the ith vertex in the ordering  */

        /* make a clique of v's neighbors in h... */
        std::vector<int> neighbors;
        for (int w : h.neighbors(v))
            neighbors.push_back(w);
        for (int j = 0; j+1 < neighbors.size(); j++) {
            int w = neighbors[j];
            for (int k = j+1; k < neighbors.size(); k++) {
                int x = neighbors[k];
                if (pi[w] > pi[v] && pi[x] > pi[v]) { /* ... if this ordering is satisfied */
                    if (!h.adjacent(w,x))
                        h.add_edge(w,x);
                }
            }
        }
    }

    return h;
}

int min_vertex(const Graph& g, std::vector<Set>& set_list, int last) {
    int min = last-1;
    while (set_list[min].empty())
        ++min;

    int u = *(set_list[min].begin());
    set_list[min].erase(u);
    return u;
}

/* O(n+n(m+(m*m/2)+m)) = O(n*m + n*m^2 + n*m) = O(n * m^2) */
std::vector<int> greedy_degree(const Graph& g, int n) {
    /* Algorithm 3 from Treewidth computations I. Upper bounds */

    Graph h = g; /* h = g */
    std::vector<Set> deg_sets(n, Set());
	int max_deg = 0;

    /* initialize deg_sets */
    for (Map<Set>::const_iterator it = h.begin(); it != h.end(); ++it) {
		int degree = it->second.size();
		deg_sets[degree].insert(it->first);
	}

    /* iterate over h */
    std::vector<int> ordering;
    for (int i = 0; i < n; i++) {
        int v = min_vertex(h, deg_sets, 1); /* choose v as a vertex of smallest degree in h */
        if (!h.contains_vertex(v)) {
	  i--;
	  continue;
	}
	ordering.push_back(v); /* let v be the ith vertex in the ordering */

        /* make a clique of v's neighbors in h */
        std::vector<int> neighbors;
        for (int w : h.neighbors(v))
            neighbors.push_back(w);
        for (int j = 0; j+1 < neighbors.size(); j++) {
            int w = neighbors[j];
            int w_deg = h.degree(w);
            int count = 0;
            for (int k = j+1; k < neighbors.size(); k++) {
                int x = neighbors[k];
                int x_deg = h.degree(x);
                if (!h.adjacent(w,x)) {
                    h.add_edge(w,x);
                    /* increase degree of x */
                    deg_sets[x_deg].erase(x);
                    deg_sets[x_deg+1].insert(x);
                    count++;
                }
            }
            /* increase degree of w */
            if (count > 0) {
                deg_sets[w_deg].erase(w);
                deg_sets[w_deg+count].insert(w);
            }
        }

        /* remove v and all its edges */
        for (int w : h.neighbors(v)) {
            int degree = h.degree(w);
    		deg_sets[degree].erase(w);
    	    deg_sets[degree-1].insert(w);
    	}
        h.remove_vertex(v);
    }

    return ordering;
}

/* O(n*m^2+(n*(m+(m*m/2*o)+m(m^2)))) = O(n*m^2 + n*m + n*m^2*o + n*m^3) = O(n*m^3) */
std::vector<int> greedy_fill_in(const Graph& g, int n) {
    /* Algorithm 3 from Treewidth computations I. Upper bounds */

    /* calculate max possible amount of fill any v could have */
    int max_fill = 0;
    for (int i = n-1; i > 0; i--)
        max_fill += i;

    Graph h = g; /* h = g */
    Map<int> fill(h.size());
    std::vector<Set> fill_sets(max_fill, Set());

    /* initialize fill and fill_sets */
    for (Map<Set>::const_iterator it = h.begin(); it != h.end(); ++it) {
		int f = fill_edges(h, it->first);
		fill_sets[f].insert(it->first);
        fill[it->first] = f;
	}

    /* iterate over h */
    std::vector<int> ordering;
    for (int i = 0; i < n; i++) {
        int v = min_vertex(h, fill_sets, 1); /* choose v as a vertex that causes the smallest number of fill edges */
        ordering.push_back(v); /* let v be the ith vertex in the ordering */

        /* make a clique of v's neighbors in h */
        std::vector<int> neighbors;
        for (int w : h.neighbors(v))
            neighbors.push_back(w);
        for (int j = 0; j+1 < neighbors.size(); j++) {
            int w = neighbors[j];
            for (int k = j+1; k < neighbors.size(); k++) {
                int x = neighbors[k];
                if (!h.adjacent(w,x)) {
                    h.add_edge(w,x);
                    /* decrease the fill of shared neighbors */
                    Set w_neighbors = h.neighbors(w);
                    Set x_neighbors = h.neighbors(x);
                    if (w_neighbors.size() <= x_neighbors.size()) {
                        for (int y : w_neighbors) {
                            if (x_neighbors.contains(y) && y != v) {
                                fill_sets[fill[y]].erase(y);
                                int new_fill = fill[y]-1;
                                fill_sets[new_fill].insert(y);
                                fill[y] = new_fill;
                            }
                        }
                    } else {
                        for (int y : x_neighbors) {
                            if (w_neighbors.contains(y) && y != v) {
                                fill_sets[fill[y]].erase(y);
                                int new_fill = fill[y]-1;
                                fill_sets[new_fill].insert(y);
                                fill[y] = new_fill;
                            }
                        }
                    }
                }
            }
        }

        /* remove v and all its edges */
        h.remove_vertex(v);
        int j = 0;
        for (int w : neighbors) {
            fill_sets[fill[w]].erase(w);
            int new_fill = fill_edges(h, w);
            fill_sets[new_fill].insert(w);
            fill[w] = new_fill;
            j++;
        }
    }

    return ordering;
}

Graph minimal_triangulation(const Graph& g, int n) {
    /* Algorithm 6 from Treewidth computations I. Upper bounds */

    Graph g_prime = g; /* g' = g */
    std::vector<Set> deg_sets(n, Set());
	int max_deg = g_prime.size()+1;

    /* initialize deg_sets */
    for (Map<Set>::const_iterator it = g_prime.begin(); it != g_prime.end(); ++it) {
		int degree = it->second.size();
		deg_sets[degree].insert(it->first);
	}

    int s = chordal(g_prime, deg_sets);
    while (s != -1) { /* while g' is not chordal */

        /* complete s in h */
        std::vector<int> neighbors;
        for (int w : g_prime.neighbors(s))
            neighbors.push_back(w);
        for (int j = 0; j+1 < neighbors.size(); j++) {
            int w = neighbors[j];
            int w_deg = g_prime.degree(w);
            int count = 0;
            for (int k = j+1; k < neighbors.size(); k++) {
                int x = neighbors[k];
                int x_deg = g_prime.degree(x);
                if (!g_prime.adjacent(w,x)) {
                    g_prime.add_edge(w,x);
                    /* increase degree of x */
                    deg_sets[x_deg].erase(x);
                    deg_sets[x_deg+1].insert(x);
                    count++;
                }
            }
            /* increase degree of w */
            if (count > 0) {
                deg_sets[w_deg].erase(w);
                deg_sets[w_deg+count].insert(w);
            }
        }
        s = chordal(g_prime, deg_sets);
    }

    return g_prime;
}

int chordal(Graph& g, std::vector<Set>& deg_sets) {
    int n = g.size();
    for (int i = 0; i < n; i++) {
        int min_v = min_vertex(g, deg_sets, 1);
        if (!clique(g, min_v))
            return min_v;

        /* remove v and all its edges */
        for (int w : g.neighbors(min_v)) {
            int degree = g.degree(w);
    		deg_sets[degree].erase(w);
    	    deg_sets[degree-1].insert(w);
    	}
        g.remove_vertex(min_v);
    }

    return -1; /* if chordal */
}

bool clique(const Graph& g, int v) {
    int n = g.neighbors(v).size();
    int v_neighbors [n];
    Set clique;

    int i = 0;
    clique.insert(v);
    for (int w : g.neighbors(v)) {
        clique.insert(w);
        v_neighbors[i] = w;
        i++;
    }

    for (int j = 0; j < n; j++) {
        int clique_size = clique.size()-1;
        int vj = v_neighbors[j];
        for (int w : g.neighbors(vj))
            if (clique.contains(w))
                clique_size--;
        if (clique_size > 0)
            return false;
    }

    return true;
}

int lowest_neighbor(const Graph& g, int v, std::vector<int> r_ordering) {
    Map<int> pi;
    for (int i = 0; i < r_ordering.size(); i++)
        pi[r_ordering[i]] = i;

    int lowest;
    for (int n : g.neighbors(v)) {
        lowest = n;
        break;
    }

    for (int n : g.neighbors(v))
        if (pi[n] > pi[lowest])
            lowest = n;

    return pi[lowest];
}

/* O(m+(m*m/2)) = O(m^2) where m = max neighbors_v */
int fill_edges(const Graph& g, int v) {
    int non_adjacent_pairs = 0;

    std::vector<int> neighbors;
    for (int w : g.neighbors(v))
        neighbors.push_back(w);

    for (int j = 0; j+1 < neighbors.size(); j++) {
        int w = neighbors[j];
        for (int k = j+1; k < neighbors.size(); k++) {
            int x = neighbors[k];
            if (!g.adjacent(w,x))
                non_adjacent_pairs++;
        }
    }

    return non_adjacent_pairs;
}

}}
