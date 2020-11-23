
#include "treewidth.hpp"

namespace sr_apx {
namespace treewidth {

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
        revdeg.resize(maxdeg + 1);
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

}}
