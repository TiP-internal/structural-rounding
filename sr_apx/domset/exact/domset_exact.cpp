
#include "sr_apx/domset/exact/domset_exact.hpp"

#include <deque>
#include <iostream>

namespace sr_apx{
namespace domset {
namespace exact {

//--- Calculating Solution
int get_soln_row_index(Table* table, const Set& optional_verts){
    // Returns the index of the row in the last table w/ smallest dominating set size.
    int soln_index=-1;
    int min_Ai = INF;
    std::vector<int> coloring;

    // finds the row w. smallest solution size.
    for(int i=0; i<table->get_table_size(); i++) {
        int curr_row_Ac = table->get_rows_Ac(i); //Ac value of row i

        if(curr_row_Ac < min_Ai) {
            int num_notdom=0;
            coloring = table->get_rowcol(i);

            for(int j=0; j<table->vertices.size(); j++) {
                int c = coloring[j];
                int xt = table->vertices[j];

                // a final coloring cant have any vertices colored as
                // not dominated (unless optional)
                if(c==NOT_DOMINATED && !optional_verts.contains(xt)) num_notdom++;
            }

            if(num_notdom==0) {
                min_Ai=curr_row_Ac;
                soln_index=i;
            }
        }
    }

    if(soln_index==-1) throw ("Solution row not found.\n");
    return soln_index;
}


int get_solution(std::vector<Table*> &tables, const Set& optional_verts) {
    // Finds the size of the final solution
    int soln_size=0;
    int tabsize = tables.size();
    Table* tab;

    for(int i=0; i<tabsize; i++) {
        tab = tables.back();
        tables.pop_back();

        int soln_index=-1;
        try {
            soln_index = get_soln_row_index(tab, optional_verts);
        } catch (const char* msg) {
            std::cerr << msg << std::endl;
            exit(0);
        }
        soln_size += tab->get_rows_Ac(soln_index);

        delete tab;
        tab=nullptr;
    }
    return soln_size;
}


int get_solution(std::vector<Table*> &tables, Set& dom_set, const Set& optional_verts) {
    /*
     * This is the second pass of the alg which constructs the soln set.
     * It iterates over tables from top to bottom.
     *
     * returns number of domset vertices added for this one compoenent
     */
    Table* parent_table = tables.back();
    tables.pop_back();

    treewidth::po_bag parent_pobag = parent_table->get_pobag();
    int soln_index=-1;
    try {
        soln_index = get_soln_row_index(parent_table, optional_verts);
    } catch (const char* msg) {
        std::cerr << msg << std::endl;
        exit(0);
    }
    int soln_size = parent_table->get_rows_Ac(soln_index);
    std::vector<int> soln_coloring = parent_table->get_rowcol(soln_index);

    // Adds the vertices which have their coloring set to IN_DOMSET in the given row.
    for(int i=0; i<parent_table->vertices.size(); i++) {
        int xt = parent_table->vertices[i];
        int c = soln_coloring[i];

        if(c==IN_DOMSET) dom_set.insert(xt);
    }

    int childl_row_ind = parent_table->get_rows_childl_table_ind(soln_index);
    int childr_row_ind = parent_table->get_rows_childr_table_ind(soln_index);

    std::deque<treewidth::po_bag> unfinished_joins;
    std::deque<int> left_inds;
    int tabsize = tables.size();  //shouldnt be more than 2

    for(int i=0; i<tabsize; i++) {
        // goal is to find the correct current row from the child table during each iteration.
        int curr_row_index;
        if(parent_pobag.num_children==2) {        //root table is a join table
            unfinished_joins.push_front(parent_pobag);
            left_inds.push_front(childl_row_ind);
        }

        Table* child_table = tables.back();
        tables.pop_back();
        treewidth::po_bag child_pobag = child_table->get_pobag();  //par_bag_index, bag_index, num_children

        if(parent_pobag.bag_index != child_pobag.parent_bag_index) {
            // if the child table's parent is not the current parent
            // then it must must be the left child of an unfinished join table
            if(unfinished_joins[0].bag_index == child_pobag.parent_bag_index) {
                parent_pobag = unfinished_joins[0];  //get most recently added
                unfinished_joins.pop_front();

                childl_row_ind = left_inds[0];
                left_inds.pop_front();
                curr_row_index = childl_row_ind;
            } else {
                throw("ERROR in getting solution\n");
            }
        } else {
            // there will be a right ind if its not a leaf
            curr_row_index = childr_row_ind;
        }
        std::vector<int> coloring = child_table->get_rowcol(curr_row_index);

        // Adds the vertices which have their coloring set to IN_DOMSET in the given row.
        for(int i=0; i<child_table->vertices.size(); i++) {
            int xt = child_table->vertices[i];
            int c = coloring[i];

            if(c==IN_DOMSET) {
                dom_set.insert(xt);
            }
        }

        childl_row_ind = child_table->get_rows_childl_table_ind(curr_row_index);
        childr_row_ind = child_table->get_rows_childr_table_ind(curr_row_index);

        delete parent_table;
        parent_table = child_table;
        parent_pobag = child_pobag;
    }

    delete parent_table;
    return soln_size;
}


int calculate(const Graph& graph, treewidth::Decomposition& decomp, Set& dom_set,
              const Set& optional_verts, Variant variant, bool construct_soln) {
    /* Dynamic programming algorithm, for dominating set on graphs
     * w/ bounded treewidth. Requires nice tree decompositions.
     * Implementation of [Alber & Niedermeier '02]
     *
     * If construct_soln is true, the solution set is constructed,
     * Calculates annotated or regular version of the dom set.
     *
     * If construct_soln is false, the tables are merged, so not all tables
     * get stored separately.
     *
     * For both annotated and regular dominating set,
     * and annotated/regular versions of independent dominating set
     * and perfect dominating set.
     */
    std::vector<Table*> tables;
    std::vector<treewidth::po_bag> postorder = decomp.get_post_order();

    // The anchors are for if we find a set of tables to
    // save--now we are saving them all (during the non-merge version).
    // Set* anchors = treedecomp_reduction(graph, decomp->components_bags, postorder);
    Set anchors;
    if(construct_soln) for(int i=0; i<decomp.components_bags.size(); i++) anchors.insert(i);

    int solutionsize=0;
    for(int i=0; i<postorder.size(); i++) {
        int bag_index = postorder[i].bag_index;
        int par_bag_index = postorder[i].parent_bag_index;

        if(decomp.components_bags[bag_index].size() != 0) {
            try{
                calculate_tables(graph, decomp.components_bags,
                                postorder, tables,
                                optional_verts,
                                anchors, variant, i);
            } catch (const char* msg) {
                std::cerr << msg << std::endl;
                exit(0);
            }
        }

        if(par_bag_index==-1 || decomp.components_bags[par_bag_index].size() == 0) {
            // empty root bag or the actual root has been reached
            // calculate solution up until this point, and discard current tables
            if(tables.size()>0) {
                if(construct_soln) {
                    try {
                        solutionsize += get_solution(tables, dom_set, optional_verts);
                    } catch (const char* msg) {
                        std::cerr << msg << std::endl;
                        exit(0);
                    }
                } else {
                    solutionsize += get_solution(tables, optional_verts);
                }
            }
        }
    }

    return solutionsize;
}

void calculate_tables(const Graph& graph, std::vector<Set> &bags,
                      std::vector<treewidth::po_bag> &postorder,
                      std::vector<Table*> &tables,
                      const Set& optional_verts, Set& anchor_tables,
                      Variant variant, int po_index) {
    // Creates the proper table type.
    treewidth::po_bag pob_current = postorder[po_index];
    Table* table = nullptr;
    bool merge;

    if(pob_current.num_children==2) {          //-----------------------JOIN bag
        int tab_ind_childr, tab_ind_childl=-1;
        treewidth::po_bag pob_rightchild, pob_leftchild;
        Table* right_child_table;
        Table* left_child_table;

        //---------------Get right child table index in tables[]
        // right child will always be the most recently added
        tab_ind_childr = tables.size()-1;
        right_child_table = tables[tab_ind_childr];
        pob_rightchild = right_child_table->get_pobag();
        if(pob_rightchild.parent_bag_index != pob_current.bag_index) {
            throw("calculate_tables: parent (right) child dont align.\n");
        }
        //----------------

        //---------------Get left child table index in tables[]
        for(int j=tables.size()-2; j>=0; j--) {
            treewidth::po_bag check = tables[j]->get_pobag();
            if(check.parent_bag_index==pob_current.bag_index) {
                tab_ind_childl = j;  //left child table index
            }
        }

        if(tab_ind_childl==-1) {
            throw("calculate_tables: parent (left) child not in tables vector.\n");
        }
        left_child_table = tables[tab_ind_childl];
        pob_leftchild = left_child_table->get_pobag();
        if(pob_leftchild.parent_bag_index != pob_current.bag_index) {
            throw("calculate_tables: parent (left) child dont align.\n");
        }
        //----------------

        if(left_child_table->get_table_size() != right_child_table->get_table_size()) {
            throw("calculate_tables: join children have different dimension\n");
        }

        if(anchor_tables.contains(pob_rightchild.bag_index)
            && anchor_tables.contains(pob_leftchild.bag_index)) {
            merge = false;
            table = join_table(right_child_table, left_child_table,
                               optional_verts, pob_current, merge);
        } else {
            // join_table(tab1, tab2, ...) -- tab2 is merged into tab1. tab1 becomes parent
            merge = true;
            table = join_table(right_child_table, left_child_table,
                               optional_verts, pob_current, merge);

            delete left_child_table;
            tables.pop_back();
            left_child_table=nullptr;
            right_child_table=nullptr;
            tables.erase(tables.begin()+tab_ind_childl);
        }
    } else if(pob_current.num_children==1) {     //-------------------either INTRODUCE or FORGET bag
        treewidth::po_bag pob_child = postorder[po_index-1];
        int child_bag_table_index = tables.size()-1;

        Set& parent_bag = bags[pob_current.bag_index];
        Set& child_bag = bags[pob_child.bag_index];

        Table* child_table = tables[child_bag_table_index];

        int v;
        if(parent_bag.size() > child_bag.size()) {        // introduce node
            // gets the introduced vertex
            for(auto it=parent_bag.begin(); it!=parent_bag.end(); it++) {
                if(!child_bag.contains(*it)) v = *it;
            }

            if(anchor_tables.contains(pob_child.bag_index)) {
                // create a new table if the child is an anchor
                merge = false;
            } else {
                // child table is now converted into the current par table
                merge = true;
                tables.pop_back();
            }
            table = intro_table(graph, child_table, child_bag,
                                optional_verts, pob_current, variant, v, merge);
        } else if(parent_bag.size() < child_bag.size()) { // forget node
            // gets the forgotten vertex
            for(auto it=child_bag.begin(); it!=child_bag.end(); it++) {
                if(!parent_bag.contains(*it)) v = *it;
            }

            if(anchor_tables.contains(pob_child.bag_index)) {
                merge = false;
            } else {
                merge = true;
                tables.pop_back();
            }
            table = forget_table(child_table, optional_verts, pob_current, variant, v, merge);
        } else {
            throw("calculate_tables: parent and child should not have same size.\n");
        }
    } else if(pob_current.num_children==0){    //----------------------------LEAF bag
        table = initialize_leaf_table(graph, bags[pob_current.bag_index],
                                      optional_verts, pob_current, variant);
    } else {
            throw("calculate_tables: wrong number of children in nice decomp.\n");
    }
    tables.push_back(table);
}


Table* init_parent_table(Table* child1_table, bool merge) {
    // how the parent table is initialized depends on whether
    // the child table is being merged into the parent or not.
    Table* par_table = nullptr;
    if(!merge) {
        // creating new table  (constructive version)
        par_table = new Table();
        par_table->vertices = child1_table->vertices; // k time
    }else if(merge) {
        // merging
        par_table = child1_table;
    }
    return par_table;
}


int locally_valid_coloring(const Graph& graph, Table* table, const Set& optional_verts,
                           int row_index, Variant variant) {
    /* Determines if the given coloring is valid or not.
     *
     * To be invalid, there must be a vertex colored as DOMINATED,
     * yet, none of its neighbors are set to IN_DOMSET.
     */
    int A_c = 0;
    const std::vector<int> row_coloring  = table->get_rowcol(row_index);
    int xt, coloring_i, num_dominators;
    bool valid=true, is_xtoptional;

    for(int i=0; i<row_coloring.size(); i++) {  //loop over the current coloring
        xt = table->vertices[i];
        coloring_i = row_coloring[i];

        if(coloring_i == IN_DOMSET) {
            if(A_c < INF) A_c++;
        } else if(coloring_i == DOMINATED) {
            is_xtoptional = optional_verts.contains(xt);
            num_dominators = 0;

            if(!is_xtoptional) {
                num_dominators = get_num_dominators(graph, table->vertices, row_coloring, xt);
                if(num_dominators<1) valid=false;
            }
            if(!valid) A_c = INF;
        }
    }
    return A_c;
}


Table* initialize_leaf_table(const Graph& graph, const Set& bag, const Set& optional_verts,
                             treewidth::po_bag pob_current, Variant variant) {
    /*
     * Creates a new leaf table.
     */
    Table* table = new Table();
    int row_update_ind, row2_ind, row3_ind;
    int tab_size, new_Ac;
    int i=0;

    for(auto it=bag.begin(); it!=bag.end(); it++) {
        int v = *it;
        table->vertices.push_back(v);

        if(i==0) {
            // initialize first three rows
            row_update_ind = table->create_row(IN_DOMSET);
            row2_ind = table->create_row(DOMINATED);
            row3_ind = table->create_row(NOT_DOMINATED);

            if(bag.size()==1) {
                // if theres only one vertex in bag, update Ac values
                tab_size = table->get_table_size();
                for(int j=0; j<tab_size; j++) {
                    new_Ac = locally_valid_coloring(graph, table, optional_verts, j, variant);
                    table->update_Ac(j, new_Ac);
                }
            }
        } else  {
            tab_size = table->get_table_size();
            for(int j=0; j<tab_size; j++) {
                row_update_ind = j;

                row2_ind = table->create_row(row_update_ind, DOMINATED);
                row3_ind = table->create_row(row_update_ind, NOT_DOMINATED);
                table->update_row_add(row_update_ind, IN_DOMSET);

                if(i==bag.size()-1) {
                    // find the Ac values once the rows are finished initializing
                    new_Ac = locally_valid_coloring(graph, table, optional_verts, row2_ind, variant);
                    table->update_Ac(row2_ind, new_Ac);

                    new_Ac = locally_valid_coloring(graph, table, optional_verts, row3_ind, variant);
                    table->update_Ac(row3_ind, new_Ac);

                    new_Ac = locally_valid_coloring(graph, table, optional_verts, row_update_ind, variant);
                    table->update_Ac(row_update_ind, new_Ac);
                }
            }
        }
        i++;
    }

    if(table->get_table_size()!=pow(3, table->vertices.size())) {
        throw("ERROR: in creating all colorings. table_size=%d \
        3^vertices.size()=%d", table->get_table_size(), pow(3, table->vertices.size()));
    }

    table->set_pobag(pob_current);
    return table;
}


int* best_join_child_rows(Table* childk, Table* childj,
                          const Set& optional_verts,
                          std::vector<int> &rowk_coloring) {
    /* In [Alber & Niedermeier '02], this function corresponds to the
     * finding the minimum c' and c'' to divide c.
     *
     * Finds the two best child rows for the current join row.
     */
    // First, find all possible c' and c'' which divide c.
    std::vector<std::vector<int>> c_prime;
    std::vector<int> c_primeprime;  //colorings for child k table

    bool rearrange=false;
    int num_ones = 0;   //number of 1's (IN_DOMSET's) in the coloring.

    // Finds if the vertices in childk and childj are in the same order
    // if they aren't, then the colorings will need to be rearranged
    for(int i=0; i<childj->vertices.size(); i++) {
        int xt = childj->vertices[i];
        int other = childk->vertices[i];
        if(xt != other) rearrange=true;
    }

    for(int i=0; i<rowk_coloring.size(); i++) {
        int c_t = rowk_coloring[i];
        int xt = childk->vertices[i];
        bool is_xtoptional = optional_verts.contains(xt);  // xt in optional set or not?

        if(c_t == IN_DOMSET) num_ones++;

        int cprime_size = c_prime.size();
        int cprimeprime_size = c_primeprime.size();

        if(cprime_size != cprimeprime_size) throw("Cannot divide colorings in minAi_c().\n");

        // is_xoptional will only be true is we want the annotated version.
        if(c_t==IN_DOMSET || c_t==NOT_DOMINATED || is_xtoptional) {
            if(cprime_size==0 && cprimeprime_size==0) {
                std::vector<int> cprime0;
                cprime0.push_back(c_t);
                c_prime.push_back(cprime0);    //c_t' == c_t

                int cprimeprime0 = c_t;
                c_primeprime.push_back(cprimeprime0);   //c_t'' == c_t
            } else {
                for(int t=0; t<cprime_size; t++) {
                    //c_t' == c_t
                    c_prime[t].push_back(c_t);

                    //c_t'' == c_t
                    c_primeprime[t] = c_primeprime[t]*10+c_t;
                }
            }
        } else if (c_t == DOMINATED ) {
            if(cprime_size==0 && cprimeprime_size==0) {
                //---- c'_t = dominated
                std::vector<int> cprime0;
                cprime0.push_back(DOMINATED);
                c_prime.push_back(cprime0);

                //c''_t = not_dominated
                int cprimeprime0 = NOT_DOMINATED;
                c_primeprime.push_back(cprimeprime0);

                //---- c'_t = not_dominated
                std::vector<int> cprime1;
                cprime1.push_back(NOT_DOMINATED);
                c_prime.push_back(cprime1);

                //c''_t = dominated
                int cprimeprime1 = DOMINATED;
                c_primeprime.push_back(cprimeprime1);
            } else {
                for(int t=0; t<cprime_size; t++) {
                    //1. make copy
                    std::vector<int> cprime_t = c_prime[t];

                    //2. update current
                    c_prime[t].push_back(DOMINATED);

                    //3. add to second pair
                    cprime_t.push_back(NOT_DOMINATED);

                    //4. add second pair to end of vector
                    c_prime.push_back(cprime_t);

                    //---------
                    //1. make copy
                    int cprimeprime_t = c_primeprime[t];

                    //2. update current
                    c_primeprime[t]=c_primeprime[t]*10+NOT_DOMINATED;

                    //3. add to second pair
                    cprimeprime_t = cprimeprime_t*10+DOMINATED;

                    //4. add second pair to end of vector
                    c_primeprime.push_back(cprimeprime_t);
                }
            }
        }
    }

    // Now, find Ai(c) ← min{Aj(c')+ Ak(c'') − #1(c) | c' and c'' divide c }
    int min_Ai = INF;
    int A_jcprime_ind_final = -1;
    int A_kprimeprime_ind_final = -1;

    for(int i=0; i<c_prime.size(); i++) {
        std::vector<int> cprime_key = c_prime[i];
        int cprimeprime_key = c_primeprime[i];

        int cprime_key_val=0;
        // the keys now need to be flipped so they align with j's coloring order
        if(rearrange) {
            //int flipped_cprime_key = flip_coloring(childj, childk->vertices, cprime_key);
            // rearranges  the coloring to match the correct childs table ordering.
            int flipped_cprime_key=-1;
            for(int i=0; i<childj->vertices.size(); i++) {  //left
                for(int k=0; k<childk->vertices.size(); k++) { //right
                    int xt = childj->vertices[i];
                    int other = childk->vertices[k];
                    if(xt == other && i==0) flipped_cprime_key=cprime_key[k];
                    else if(xt==other) flipped_cprime_key = flipped_cprime_key*10+cprime_key[k];
                }
            }
            cprime_key_val = flipped_cprime_key;
        } else {
            int key=-1;
            for(int i=0; i<childj->vertices.size(); i++) {
                int xt = childj->vertices[i];
                if(i==0) key=cprime_key[i];
                else key = key*10+cprime_key[i];
            }
            cprime_key_val=key;
        }

        int A_jcprime_ind = childj->lookup_table_index(cprime_key_val);
        int A_jcprime = childj->lookup_Ac(cprime_key_val);

        int A_kprimeprime_ind = childk->lookup_table_index(cprimeprime_key);
        int A_kprimeprime = childk->lookup_Ac(cprimeprime_key);

        int val;
        if(A_jcprime==INF || A_kprimeprime==INF) val = INF;
        else val = A_jcprime + A_kprimeprime - num_ones;

        if(val <= min_Ai) {
            min_Ai = val;
            A_jcprime_ind_final = A_jcprime_ind;
            A_kprimeprime_ind_final = A_kprimeprime_ind;
        }
    }

    int *min_vals = new int[3];
    min_vals[0] = min_Ai;
    min_vals[1] = A_jcprime_ind_final;      // j is left child table
    min_vals[2] = A_kprimeprime_ind_final;  // k is right

    return min_vals;
}


Table* join_table(Table* child1_table, Table* child2_table,
                  const Set& optional_verts, treewidth::po_bag pob_current, bool merge) {
    /* If merge is true, the par_table is set to the child_table
     * If merge if false, set par_table = new Table();,
     */
    Table* par_table = init_parent_table(child1_table, merge);
    int tab_size = child2_table->get_table_size();
    int rowk_index;
    std::vector<int> rowk_coloring;
    int* min_vals;

    for(int i=0; i<tab_size; i++) {  //4^k
        if(!merge) {  //creating new table  (constructive version)
            rowk_index = par_table->copyin_row(child1_table, i);
        } else rowk_index = i;

        rowk_coloring = par_table->get_rowcol(rowk_index);

        try {
            min_vals = best_join_child_rows(child1_table, child2_table, optional_verts, rowk_coloring);
        } catch (const char* msg) {
            std::cerr << msg << std::endl;
            exit(0);
        }

        par_table->update_Ac(rowk_index, min_vals[0]);
        par_table->update_rows_childl_table_ind(rowk_index, min_vals[1]);
        par_table->update_rows_childr_table_ind(rowk_index, min_vals[2]);

        delete[] min_vals;
    }
    par_table->set_pobag(pob_current);
    return par_table;
}


Table* intro_table(const Graph& graph, Table* child_table, const Set& child_bag, const Set& optional_verts,
                   treewidth::po_bag pob_current, Variant variant, int v, bool merge) {
    /*
     * Creates a table that includes the introduced vertex.
     */
    Table* parent_table = init_parent_table(child_table, merge);
    parent_table->vertices.push_back(v);
    Set neighbors_v;
    for (int x : graph.neighbors(v)) {
        if (child_bag.contains(x)) {
            neighbors_v.insert(x);
        }
    }

    bool is_xoptional = optional_verts.contains(v);
    int table_size = child_table->get_table_size();

    int row1_ind, row2_ind, row3_ind;
    int row1_key, row2_key, row3_key;

    int A_phi, best_child_key;
    int current_row_ind=-1;

    for(int i=0; i<table_size; i++) {  // looping over child table
        // initialize the new introduced rows.
        row1_ind = parent_table->copyin_introrow(child_table, i, IN_DOMSET);
        row1_key = parent_table->get_rows_key(row1_ind);

        row2_ind = parent_table->copyin_introrow(child_table, i, DOMINATED);
        row2_key = parent_table->get_rows_key(row2_ind);

        row3_ind = parent_table->copyin_introrow(child_table, i, NOT_DOMINATED);
        row3_key = parent_table->get_rows_key(row3_ind);

        if(merge) {
            current_row_ind++;
            parent_table->update_table_lookups(row1_key, current_row_ind);

            current_row_ind++;
            parent_table->update_table_lookups(row2_key, current_row_ind);

            current_row_ind++;
            parent_table->update_table_lookups(row3_key, current_row_ind);
        }

        int num_dominators = 0;
        const std::vector<int> rupdate_coloring = parent_table->get_rowcol(row1_ind);
        const std::vector<int> row2_coloring = parent_table->get_rowcol(row2_ind);

        // -----Introducing a vertex w/ coloring set to IN_DOMSET.
        // r_update: x is IN_DOMSET=1
        best_child_key = best_intro_child_row(neighbors_v,
                                              child_table->vertices,
                                              rupdate_coloring, v);   //k--phi()

        // domset variant
        A_phi = child_table->lookup_Ac(best_child_key);
        if(A_phi != INF) A_phi++;

        // Update r_update row
        int index = child_table->lookup_table_index(best_child_key);
        parent_table->update_rows_childl_table_ind(row1_ind, -1);
        parent_table->update_rows_childr_table_ind(row1_ind, index);
        parent_table->update_Ac(row1_ind, A_phi);


        //------Introducing a vertex w/ coloring set to DOMINATED.
        // r2: x is DOMINATED=0
        if(!is_xoptional) {
            // check that the introduced vertex w. the coloring of DOMINATED is justified.
            for(int j=0; j<row2_coloring.size(); j++) {
                int x = parent_table->vertices[j];

                // x is set to DOMINATED and has a neighbor IN_DOMSET
                if(row2_coloring[j]==IN_DOMSET
                    && neighbors_v.contains(x)) num_dominators++;
            }
            if(num_dominators<1) {  // not valid coloring
                parent_table->update_Ac(row2_ind, INF);
            }
        }

        //------Introducing a vertex w/ coloring set to NOT_DOMINATED.
        // r3: x is NOT_DOMINATED=3  -- do nothing
    }

    if(merge) {
        // it's appending all rows to the end of parent. so if we are merging,
        // we will need to delete the previous rows from the child
        for(int i=0; i<table_size; i++) {
            int key = parent_table->get_rows_key(0);
            parent_table->pop_front_row();
            parent_table->table_lookups_remove(key);
        }
    }

    parent_table->set_pobag(pob_current);
    return parent_table;
}


Table* forget_table(Table* child_table, const Set& optional_verts, treewidth::po_bag pob_current,
                    Variant variant, int v, bool merge) {
    /*
     * Creates table with the forgotten vertex v.
     */
    Table* par_table = init_parent_table(child_table, merge);
    int v_index = par_table->get_vertex_col_index(v);  // k
    int table_size = child_table->get_table_size();
    bool is_xoptional = optional_verts.contains(v);  // v in optional set or not?

    par_table->vertices.erase(par_table->vertices.begin()+v_index); // delete v from verts vec. k

    for(int j=0; j<table_size; j++) {   //  loop over child table    //3^k
        int row_ind_child=j;     // current child row index
        int true_row_ind_child=j;  // tracks the actual row index in child table
        int row_key_child=child_table->get_rows_key(row_ind_child);

        //---- creates the key for if the child was a par row
        int row_key_par=-99, forgotten_color;
        std::vector<int> coloring_child = child_table->get_rowcol(row_ind_child);
        std::vector<int> coloring_par;
        for(int i=0; i<coloring_child.size(); i++) {
            if(row_key_par==-99 && i!=v_index) row_key_par=coloring_child[i];
            else if(i!=v_index) {
                row_key_par=row_key_par*10+coloring_child[i];
            }
            if(i!=v_index) coloring_par.push_back(coloring_child[i]);
            if(i==v_index) forgotten_color=coloring_child[i];
        }

        //---- get row_ind_par
        int row_ind_par=par_table->lookup_table_index(row_key_par);

        // if the current child row not in parent
        if(row_ind_par==-1) {
            row_ind_par = par_table->copyin_forgetrow(child_table, row_ind_child, v_index);
            par_table->update_rows_childl_table_ind(row_ind_par, -1);
            par_table->update_rows_childr_table_ind(row_ind_par, true_row_ind_child);
        }
        //-----

        // now potentially update the current parent row
        int Ac_par, Ac_child;
        Ac_par = par_table->get_rows_Ac(row_ind_par);
        Ac_child = child_table->get_rows_Ac(row_ind_child);

        if(forgotten_color!=NOT_DOMINATED && Ac_child < Ac_par) {
            par_table->update_rows_childl_table_ind(row_ind_par, -1);
            par_table->update_rows_childr_table_ind(row_ind_par, true_row_ind_child);
            par_table->update_Ac(row_ind_par, Ac_child);
        }
    }

    if(merge) { // NOTE  this is probably not ideal
        // it's appending all rows to the end of parent. so if we are merging,
        // we will need to delete the previous rows from the child
        for(int i=0; i<table_size; i++) {
            int key = par_table->get_rows_key(0);
            par_table->pop_front_row();
            par_table->table_lookups_remove(key);
        }

        for(int i=0; i<par_table->get_table_size(); i++) {
            int curr_key = par_table->get_rows_key(i);
            par_table->update_table_lookups(curr_key, i);
        }
    }

    par_table->set_pobag(pob_current);
    return par_table;
}


//----- Helper functions
int best_intro_child_row(const Set& neighbors, std::vector<int> &vertices,
                         const std::vector<int> &coloring, int introduced_v) {
   /* In [Alber & Niedermeier '02], this corresponds to phi(c) in introduce nodes.
    *
    *  φ : {0, ˆ0, 1}^nj → {0, ˆ0, 1}^nj
    * on the set of colorings of Xj.
    * For c =(c1,... ,c_nj ) ∈ {0, ˆ0, 1}^nj , let φ(c):= (c'_1,... ,c'_nj )
    * such that
    *
    * c'_t = 0ˆ if t ∈ {p1,... ,ps} and ct = 0  OR
    * c'_t = c_t otherwise.
    *
    * returns key of the best child coloring.
    */
    int new_col;
    for(int i=0; i<vertices.size(); i++) {
        int curr_v = vertices[i];
        if(curr_v != introduced_v) {
            int curr_v_c = coloring[i];

            // the introduce vertex has a neighbor w. color set to DOMINATED.
            if(curr_v_c==DOMINATED && neighbors.contains(curr_v)) {
                if(i==0) new_col = NOT_DOMINATED;
                else new_col = new_col*10+NOT_DOMINATED;
            } else {
                if(i==0) new_col = curr_v_c;
                else new_col = new_col*10+curr_v_c;
            }
        }
    }
    return new_col;
}

int get_num_dominators(const Graph& graph, std::vector<int> &vertices,
                       const std::vector<int> &coloring, int target_vert) {
    //returns number of vertices dominating target_vert.
    int num_dominators=0;

    for(int j=0; j<vertices.size(); j++) {  //at most k
        int coloring_j = coloring[j];

        if(vertices[j]!=target_vert && coloring_j== IN_DOMSET) {
            if(graph.adjacent(vertices[j], target_vert)) {
                num_dominators++;
            }
        }
    }
    return num_dominators;
}

}}}
