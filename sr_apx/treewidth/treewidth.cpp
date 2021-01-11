
#include "treewidth.hpp"

namespace sr_apx::treewidth {

Decomposition::Decomposition() {
    tw = 0;
}

Decomposition::Decomposition(const Graph& graph): Decomposition() {
    build_decomposition(graph);
}

int Decomposition::treewidth() {
    return tw - 1;
}

std::vector<po_bag> Decomposition::get_post_order() {
    std::vector<po_bag> po;
    for (auto it = pre_order.rbegin(); it != pre_order.rend(); ++it) {
        po.push_back(*it);
    }

    return po;
}

// returns placement of new bag
// parent_index should be 0 for root bag
// only modifies pre_order and components_bags
int Decomposition::add_bag(int parent, bool last_child, Set& bag) {
    // ensures that root bag exists
    if (components_bags.empty()) {
        components_bags.push_back(Set());
        pre_order.push_back(po_bag(0, -1));
    }

    if (bag == components_bags[parent]) {
        printf("%s\n", "equal bags");
        return parent;
    }

    int true_parent = pre_order[parent].current_join_child;

    if (pre_order[true_parent].num_children > 0 && !last_child) {
        // create a new join bag to hold the child
        int index = components_bags.size();
        components_bags.push_back(components_bags[true_parent]);
        pre_order.push_back(po_bag(index, true_parent));
        pre_order[true_parent].num_children++;
        pre_order[parent].current_join_child = index;
        true_parent = index;
    }

    if (pre_order[true_parent].num_children > 0 || !last_child) {
        int index = components_bags.size();
        components_bags.push_back(components_bags[true_parent]);
        pre_order.push_back(po_bag(index, true_parent));
        pre_order[true_parent].num_children++;
        true_parent = index;
    }

    // create the child
    int use_parent = true_parent;
    Set use_bag = components_bags[true_parent];

    // add introduce bags
    for (int x : components_bags[true_parent]) {
        // only look at vertices to be added
        if (bag.contains(x)) {
            continue;
        }

        use_bag.erase(x);
        int index = components_bags.size();
        components_bags.push_back(use_bag);
        pre_order.push_back(po_bag(index, use_parent));
        pre_order[use_parent].num_children++;
        use_parent = index;
    }

    // add forget bags
    for (int x : bag) {
        // only look at vertices to be removed
        if (components_bags[true_parent].contains(x)) {
            continue;
        }

        use_bag.insert(x);
        int index = components_bags.size();
        components_bags.push_back(use_bag);
        pre_order.push_back(po_bag(index, use_parent));
        pre_order[use_parent].num_children++;
        use_parent = index;
    }

    tw = tw < bag.size() ? bag.size() : tw;

    // return index of final bag
    return components_bags.size() - 1;
}

void Decomposition::build_decomposition(const Graph& graph) {
    std::vector<Set> components = graph.connected_components();
    int nc = components.size();

    for (int i = 0; i + 1 < nc; ++i) {
        Graph component = graph.subgraph(components[i]);
        Set empty;
        tree_decomp(component, components[i], empty, 0, false);
    }

    Graph component = graph.subgraph(components[nc - 1]);
    Set empty;
    tree_decomp(component, components[nc - 1], empty, 0, true);
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
    for (int w : g.neighbors(vn))
        g.remove_edge(vn,w);
    g.remove_vertex(vn);
    bags[n-1] = x_vn;
    add_bag(0, false, x_vn);

    /* add the rest of the bags */
    for (int i = n-2; i > -1; i--) {
        Set x_vi;
        int vi = ordering[i];
        x_vi.insert(vi); /* add vi to bag */

        /* add vi's forward neighbors */
        for (int w : g.neighbors(vi)) {
            x_vi.insert(w);
            g.remove_edge(vi,w);
        }
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

    if (Z_union_W.size() <= tw) {
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

// end of Decomposition object //////////////////////////

Set vertex_delete(const Graph& graph, int w) {
    /*
     * TreeWidthNodeEdit() from SR paper.
     *
     */

    Decomposition decomp(graph);
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
    std::vector<Set> components = sub_g.connected_components();

    for (int i = 0; i < components.size(); ++i) {
        Graph component = graph.subgraph(components[i]);
        Set edited = vertex_delete(component, w);
        S.insert(edited.begin(), edited.end());
    }

    return S;
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

    while (A_count + C_count < W.size() - ((2 * W.size()) / 3) || !revdeg[0].empty() || (!revdeg[1].empty() && A_count + C_count < W.size() - 1)) {
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


// greedy heuristics

int chordal(const Graph&, std::vector<Set>&);
bool clique(const Graph&, int);
int lowest_neighbor(const Graph&, int, std::vector<int>);
int min_vertex(const Graph&, std::vector<Set>&);
int fill_edges(const Graph&, int);
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
        int v = ordering[i]; /* let v be the ith vertex in the ordering*/

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
    std::vector<Set> deg_sets;
	int max_deg = 0;

    /* initialize deg_sets */
    for (Map<Set>::const_iterator it = h.begin(); it != h.end(); ++it) {
		int degree = it->second.size();
		max_deg = degree > max_deg ? degree : max_deg;
		deg_sets.resize(max_deg+1);
		deg_sets[degree].insert(it->first);
	}

    deg_sets.resize(max_deg * 1000); // TODO

    /* iterate over h */
    int last_deg = 1;
    std::vector<int> ordering;
    for (int i = 0; i < n; i++) {
        int v = min_vertex(h, deg_sets, last_deg); /* choose v as a vertex of smallest degree in h */
        ordering.push_back(v); /* let v be the ith vertex in the ordering */
        last_deg = h.degree(v); /* new min degree vertex can be no less than last_deg-1 */

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
                    max_deg = x_deg+1 > max_deg ? x_deg+1 : max_deg;
                    deg_sets.resize(max_deg+1);
                    deg_sets[x_deg+1].insert(x);
                    count++;
                }
            }
            /* increase degree of w */
            if (count > 0) {
                deg_sets[w_deg].erase(w);
                max_deg = w_deg+count > max_deg ? w_deg+count : max_deg;
                deg_sets.resize(max_deg+1);
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

    Graph h = g; /* h = g */
    Map<int> fill(h.size());
    std::vector<Set> fill_sets;
	int max_fill = 0;

    /* initialize fill and fill_sets */
    for (Map<Set>::const_iterator it = h.begin(); it != h.end(); ++it) {
		int f = fill_edges(h, it->first);
		max_fill = f > max_fill ? f : max_fill;
		fill_sets.resize(max_fill+1);
		fill_sets[f].insert(it->first);
        fill[it->first] = f;
	}

    fill_sets.resize(max_fill*1000); // TODO

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
            max_fill = new_fill > max_fill ? new_fill : max_fill;
            fill_sets.resize(max_fill+1);
            fill_sets[new_fill].insert(w);
            fill[w] = new_fill;
            j++;
        }
    }

    return ordering;
}

Graph minimal_triangulation(const Graph& g) {
    /* Algorithm 6 from Treewidth computations I. Upper bounds */

    Graph g_prime = g; /* g' = g */
    std::vector<Set> deg_sets;
	int max_deg = 0;

    /* initialize deg_sets */
    for (Map<Set>::const_iterator it = g_prime.begin(); it != g_prime.end(); ++it) {
		int degree = it->second.size();
		max_deg = degree > max_deg ? degree : max_deg;
		deg_sets.resize(max_deg+1);
		deg_sets[degree].insert(it->first);
	}

    deg_sets.resize(max_deg*10000); // TODO

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
                    max_deg = x_deg+1 > max_deg ? x_deg+1 : max_deg;
                    deg_sets.resize(max_deg+1);
                    deg_sets[x_deg+1].insert(x);
                    count++;
                }
            }
            /* increase degree of w */
            if (count > 0) {
                deg_sets[w_deg].erase(w);
                max_deg = w_deg+count > max_deg ? w_deg+count : max_deg;
                deg_sets.resize(max_deg+1);
                deg_sets[w_deg+count].insert(w);
            }
        }

        s = chordal(g_prime, deg_sets);
    }

    return g_prime;
}

int chordal(const Graph& g, std::vector<Set>& deg_sets) {
    Graph h = g;
    int n = h.size();

    int last_deg = 1;
    for (int i = 0; i < n; i++) {
        int min_v = min_vertex(h, deg_sets, last_deg);
        last_deg = h.degree(min_v);

        if (!clique(h, min_v))
            return min_v;

        /* remove v and all its edges */
        for (int w : h.neighbors(min_v)) {
            int degree = h.degree(w);
    		deg_sets[degree].erase(w);
    	    deg_sets[degree-1].insert(w);
    	}
        h.remove_vertex(min_v);
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
}
