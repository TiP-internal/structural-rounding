
#include "domset_exact.hpp"

#include <cstdio>


//--- Calculating Solution
int get_soln_row_index(Table* table){
    //Returns the index of the row in the last table w/ smallest dominating set size.
    int soln_index=-1;
    int min_Ai = INF;

    //finds the row w. smallest solution size.
    for(int i=0; i<table->get_table_size(); i++) {
        Row* row = table->get_row(i);

        bool valid=true;
        for(int j=0; j<row->coloring.size(); j++) {
            int c=row->coloring[j];
            if(c==NOT_DOMINATED) valid=false;
        }

        if(valid && row->get_Ac() < min_Ai) {
            min_Ai=row->get_Ac();
            soln_index=i;
        }
    }

    return soln_index;
}


void add_to_solution(Set* dom_set, Row* row, std::vector<int> &vertices) {
    //Adds the vertices which have their coloring set to IN_DOMSET in the given row.
    for(int i=0; i<vertices.size(); i++) {
        int xt = vertices[i];
        int c = row->coloring[i];

        if(c==IN_DOMSET) dom_set->insert(xt);
    }
}


int get_solution(Table* table) {
    /*
     * Finds the size of the final solution for the optimization version of DOMSet
     */
    int soln_index = get_soln_row_index(table);

    if(soln_index==-1) printf("ERROR: soln row not found.\n");

    return table->get_row(soln_index)->get_Ac();
}


void get_solution(std::vector<Table*> &tables, Set* dom_set) {
    /*
     * This is the second pass of the alg which constructs the soln set.
     * It iterates over tables from top to bottom.
     */
    Table* final_table = tables.back();
    tables.pop_back();

    int soln_index = get_soln_row_index(final_table);
    Row* soln_row = final_table->get_row(soln_index);
    int soln_size = soln_row->get_Ac();

    add_to_solution(dom_set, soln_row, final_table->vertices);

    int childl_ind = soln_row->get_childl_table_ind();
    int childr_ind = soln_row->get_childr_table_ind();

    delete final_table;

    std::vector<int> rightinds;

    for(int i=tables.size()-1; i>=0; i--) {
        Row* curr_row;
        Table* curr_table = tables.back();
        tables.pop_back();

        if(childl_ind!=-1 && childr_ind!=-1) {          //join table
            rightinds.push_back(childr_ind);
            curr_row = curr_table->get_row(childl_ind);
        } else if(childl_ind==-1 && childr_ind==-1) {   //leaf table
            childr_ind = rightinds.back();
            rightinds.pop_back();
            curr_row = curr_table->get_row(childr_ind);
        } else if(childl_ind!=-1 && childr_ind==-1) {   //intro or forget
            curr_row = curr_table->get_row(childl_ind);
        }

        if(dom_set->size() < soln_size) add_to_solution(dom_set, curr_row, curr_table->vertices);

        childl_ind = curr_row->get_childl_table_ind();
        childr_ind = curr_row->get_childr_table_ind();

        delete curr_table;
    }
}


Set* construct_domset(Graph* graph, TreeDecomp* decomp,
                      Set* optional_verts, Variant variant) {
    /*
     * Calculates the minimum dominating set given a tree decomp.
     * This version constructs the solution set.
     *
     * Must calculate for each component in the decomp.
     * Calculates annotated or regular version of the dom set.
    */
    std::vector<std::vector<po_bag>> postorder = decomp->get_post_order();

    Set* dom_set = new Set();

    //loop over components
    for(int j=0; j<decomp->components_bags.size(); j++) {
        std::vector<Table*> tables;

        Set* anchors = treedecomp_reduction(graph, decomp->components_bags[j], postorder[j]);

        calculate_tables(graph, decomp->components_bags[j],
                         postorder[j], tables, optional_verts,
                         anchors, variant);

        get_solution(tables, dom_set);

        delete anchors;
    }
    return dom_set;
}


int calc_min_domset(Graph* graph, TreeDecomp* decomp,
                    Set* optional_verts, Variant variant) {
    /*
     * Used for the optimization version. Finds the size of the smallest dom set.
     */

    std::vector<std::vector<po_bag>> postorder = decomp->get_post_order();
    int solnsize=0;

    //loop over components
    for(int j=0; j<decomp->components_bags.size(); j++) {
        Set* empty=new Set();
        std::vector<Table*> tables;

        calculate_tables(graph, decomp->components_bags[j],
                         postorder[j], tables, optional_verts,
                         empty, variant);
        delete empty;

        solnsize += get_solution(tables[tables.size()-1]); //soln size for each component
    }
    return solnsize;
}


//---- Constructive Version
Set* treedecomp_reduction(Graph* graph, std::vector<Set*> &bags,
                          std::vector<po_bag> postorder) {
    /*
     * Tree Decomposition reduction rules to reduce the number of stored tables.
     *
     * Returns a set of integers, corresponding to bag indices which are the 'anchor' bags,
     * meaning the tables created for these bags must be stored. All other tables can be merged
     * (ie. not stored separately).
     *
     * Method from [Betzler, Neidermeier, & Uhlmann, '04].
     *
     * Not passing postorder as reference because we need to modify it here.
     */
    Set* reduced_bags = new Set();

    Set* removed_elems = new Set();
    Set* removed_nodes = new Set();

    bool cont = true;
    while(cont) {
        bool rule1=true, rule2=true, rule3=true, rule4=true, rule5=true;

        /* Rule 1: if bag Xi contains at least one element that is not in any other bag,
         * and the element is not contained in the removed elemens set,
         * put Xi in solution and add its elements to the removed elemens set.
         */
        rule1 = is_exclusive_to_singlebag(bags, reduced_bags, removed_elems, removed_nodes);


        /* Rule 2: If there is an edge e = {i, j} ∈ F such that
         *          (1) Xi ⊆ Xj and
         *          (2) w(Xi) ≥ w(Xj)
         * then connect each neighboring node of i (except for j) with j
         * and delete i together with its incident edges from T.
         */
        rule2 = is_child_subset(bags, postorder, removed_elems, removed_nodes);


        /* Rule 3: If for any bag Xi the fact a ∈ Xi implies b ∈ Xi, then
         * remove b from the bags of the tree decomposition.
         */
        rule3 = is_a_implies_b(graph, bags, reduced_bags, removed_elems, removed_nodes);


        /* Rule 4: If there is an edge e = {i, j} ∈ F such that Xi ∩ Xj = ∅,
         * then remove this edge from F.
         */
        rule4 = is_empty_bag_intersect(bags, postorder, removed_elems);


        /* Rule 5: If for a bag Xi there is a set of nodes
         * j1, j2, . . . , jl ∈ I, such that
         *      (1) Xi ⊆ UNION k=1 to l of X_jk,
         *      (2) w(Xi) ≥  SUM k=1 to l of w(X_jk),
         *      (3) the neighbors N(i) of i form a tree-like subset collection,
         *
         * then remove i from I and connect its neighbors such
         * that the tree-likeness is not violated
         */
        rule5 = is_treelike_subcollection(bags, postorder, removed_elems, removed_nodes);

        //if no rules reduced, quit.
        if(!rule1 && !rule2 && !rule3 && !rule4 && !rule5) cont=false;
    }

    delete removed_elems, removed_nodes;
    return reduced_bags;
}


void calculate_tables(Graph* graph, std::vector<Set*> &bags,
                      std::vector<po_bag> &postorder, std::vector<Table*> &tables,
                      Set* optional_verts, Set* anchor_tables, Variant variant) {
    /*
     * Dynamic programming algorithm, for dominating set on graphs
     * w/ bounded treewidth.
     *
     * Constructive version and does not reuse (all) tables. To save memory, it
     * will eventually only store anchor tables.
     *
     * For NICE tree decompositions. For both annotated and regular dominating set.
     */
    std::vector<int> table_bag_indices;  //the bag indices of the current bags in the table.
    for(int i=0; i<postorder.size(); i++) { // O(n)
        int bag_index = postorder[i].bag_index;
        int num_children = postorder[i].num_children;
        int parent_bag_index = postorder[i].parent_bag_index;

        if(num_children==2) {          //-----------------------JOIN bag
            int ind_childr = tables.size()-1;  //right child will always be the most recently added
            int ind_childl = -99;

            //NOTE needs testing
            for(int j=table_bag_indices.size()-2; j>=0; j--) {
                int ind=table_bag_indices[j];
                if(!anchor_tables->contains(ind) && ind_childl == -99) {
                    ind_childl=j;
                    continue;
                }
            }

            int bag_ind_childl = table_bag_indices[ind_childl];
            int bag_ind_childr = table_bag_indices[ind_childr];

            Table* left_child_table = tables[ind_childl];
            Table* right_child_table = tables[ind_childr];

            //update table here
            if(anchor_tables->contains(bag_ind_childr) && anchor_tables->contains(bag_ind_childl)) {
                //NOTE needs testing
                //if both children are anchor tables, create new table and add to tables
                Table* table = update_join_table(left_child_table, right_child_table, optional_verts, variant);
                tables.push_back(table);
            } else if(anchor_tables->contains(bag_ind_childr)){
                //NOTE needs testing
                //if right child is an anchor, merge left child into parent join table, dont delete
                merge_join_table(left_child_table, right_child_table,
                                 optional_verts, variant);           // reuses the right child table
                table_bag_indices.pop_back();
            } else if(anchor_tables->contains(bag_ind_childl)) {
                //NOTE needs testing
                //if left child is an anchor, merge right child into parent join, dont delete
                merge_join_table(right_child_table, left_child_table,
                                 optional_verts, variant);           // reuses the right child table
                table_bag_indices.erase(table_bag_indices.begin()+ind_childl);
            } else {
                //if neither are anchors, merge right child into parent table, delete left
                merge_join_table(right_child_table, left_child_table,
                                 optional_verts, variant);           // reuses the right child table
                tables.erase(tables.begin()+ind_childl);            // deletes the left child table.

                table_bag_indices.erase(table_bag_indices.begin()+ind_childl);  //left
                table_bag_indices.pop_back(); //right
            }
            table_bag_indices.push_back(bag_index);

        } else if(num_children==1) {     //-------------------either INTRODUCE or FORGET bag
            int child_bag_index = postorder[i-1].bag_index;
            int child_bag_table_index = tables.size()-1;

            Set* parent_bag = bags[bag_index];
            Set* child_bag = bags[child_bag_index];

            string type;
            int v;

            if(parent_bag->size() > child_bag->size()) {        // introduce node
                for(auto it=parent_bag->begin(); it!=parent_bag->end(); it++) {
                    int u = *it;
                    if(!child_bag->contains(u)) v = u;
                }

                if(anchor_tables->contains(child_bag_index)) {
                    //create a new table if the child is an anchor
                    Table* table = update_introduce_table(graph, tables[child_bag_table_index],
                                        child_bag, optional_verts, v, variant);
                    tables.push_back(table);
                } else {
                    Table* child_table = tables[child_bag_table_index];
                    merge_introduce_table(graph, child_table, child_bag, optional_verts, v, variant);
                    table_bag_indices.pop_back();
                }
                table_bag_indices.push_back(bag_index);

            } else if(parent_bag->size() < child_bag->size()) { // forget node

                for(auto it=child_bag->begin(); it!=child_bag->end(); it++) {
                    int u = *it;
                    if(!parent_bag->contains(u)) v = u;
                }

                if(anchor_tables->contains(child_bag_index)) {
                    //create a new table if the child is an anchor
                    Table* table = update_forget_table(tables[child_bag_table_index],
                                                        optional_verts, v, variant);
                    tables.push_back(table);
                } else {
                    Table* child_table = tables[child_bag_table_index];
                    merge_forget_table(child_table, optional_verts, v, variant);
                    table_bag_indices.pop_back();
                }
                table_bag_indices.push_back(bag_index);
            } else {
                printf("ERROR: in find_bagtype(), parent and child should not have same size.\n");
            }
        } else if(num_children==0){    //----------------------------LEAF bag
            //leaf bag, need to initialize table.
            //create new leaf for both anchors and non-anchors
            Table* table = initialize_leaf_table(graph, bags[bag_index], optional_verts, variant);
            tables.push_back(table);
            table_bag_indices.push_back(bag_index);
        } else {
                printf("ERROR in number of children in nice decomp.\n");
        }
    }

    //print_tables(tables);
}


Table* initialize_leaf_table(Graph* graph, Set* bag, Set* optional_verts, Variant variant) {
    /*
     * k*3^k
     *
     * Leaf tables can have more than 1 vertex.
     */
    Table* table = new Table();

    Row* r1;
    Row* r2;
    Row* r3;

    //intitial rows for first vert.
    if(bag->size()>0) {
        r1 = table->create_row(IN_DOMSET);
        r2 = table->create_row(DOMINATED);
        r3 = table->create_row(NOT_DOMINATED);
    }

    int i=0;
    for(auto it=bag->begin(); it!=bag->end(); it++) {  //3^k*k time
        int v = *it;
        table->vertices.push_back(v);

        if(i==0 && bag->size()==1) {
            //only one vertex in bag, need to find if valid.
            int tab_size = table->get_table_size();
            for(int j=0; j<tab_size; j++) {  //at most 3^ni
                Row* r_update = table->get_row(j);
                r_update->update_Ac(locally_valid_coloring(graph, optional_verts,
                                    r_update, table->vertices, variant));   //k^2
            }
        } else if(i>0) { //first vert col already initialized in constructor
            int tab_size = table->get_table_size();
            for(int j=0; j<tab_size; j++) {        //3^k
                Row* r_update = table->get_row(j);

                Row* r2 = table->create_row(r_update, DOMINATED);
                if(r2->coloring.size()==bag->size()) {
                    r2->update_Ac(locally_valid_coloring(graph, optional_verts,
                                  r2, table->vertices, variant));               //k^2
                }

                Row* r3 = table->create_row(r_update, NOT_DOMINATED);
                if(r3->coloring.size()==bag->size()) {
                    r3->update_Ac(locally_valid_coloring(graph, optional_verts,
                                  r3, table->vertices, variant));               //k^2
                }

                table->update_row_add(r_update, IN_DOMSET);
                if(r_update->coloring.size()==bag->size()) {
                    r_update->update_Ac(locally_valid_coloring(graph, optional_verts,
                                        r_update, table->vertices, variant));   //k^2
                }
            }
        }
        i++;
    }

    if(table->get_table_size()!=pow(3, table->vertices.size())) {
        printf("ERROR: in creating all colorings. table_size=%d \
        3^vertices.size()=%d", table->get_table_size(), pow(3, table->vertices.size()));
    }
    return table;
}


Table* update_introduce_table(Graph* graph, Table* child_table,
                                             Set* child_bag, Set* optional_verts,
                                             int v, Variant variant) {
    /*
     * Initializes the parent_table. Constructive version
     */
    Table* parent_table = new Table();

    parent_table->vertices = child_table->vertices; //k
    parent_table->vertices.push_back(v);

    Set* neighbors_v = child_bag->set_intersection(graph->neighbors(v));
    int table_size = child_table->get_table_size();

    for(int i=0; i<table_size; i++) {           //3^ni
        Row* r_update = new Row(child_table->get_row(i));
        parent_table->insert_row(r_update);
        r_update->set_childl_table_ind(i);

        Row* r2 = parent_table->create_row(r_update, DOMINATED);
        Row* r3 = parent_table->create_row(r_update, NOT_DOMINATED);

        //r_update: x is IN_DOMSET=1
        parent_table->update_row_add(r_update, IN_DOMSET);
        intro_vert_indomset_update(graph, child_table, neighbors_v, r_update, v, variant);

        //r2: x is DOMINATED=0
        intro_vert_dominated_update(graph, child_table, neighbors_v, optional_verts,
                                    r2, r_update, v, variant);

        //r3: x is NOT_DOMINATED=3  -- do nothing
    }
    delete neighbors_v;

    return parent_table;
}


Table* update_forget_table(Table* child_table, Set* optional_verts, int v, Variant variant) {
    /*
     * Initializes parent table.
     */
    Table* parent_table = new Table();

    int v_index = child_table->get_vertex_col_index(v); //k
    int table_size = child_table->get_table_size();
    bool is_xoptional = optional_verts->contains(v); // in optional set or not?

    parent_table->vertices = child_table->vertices; //k  NOTE could be updated
    //delete v from verts vec. k
    parent_table->vertices.erase(parent_table->vertices.begin()+v_index);

    for(int j=0; j<table_size; j++) {                //3^ni
        Row* rj = new Row(child_table->get_row(j));  //k

        int color_v = rj->coloring[v_index];
        int Aj = rj->get_Ac();
        int rj_tabind = child_table->lookup_table_index(rj->get_key());

        rj->remove_from_coloring(v_index);           //k

        if(parent_table->lookup_table_index(rj->get_key()) == -1) { //not inserted yet
            parent_table->insert_row(rj);
        }

        Row* ri = parent_table->lookup_row(rj->get_key());
        int Ai = parent_table->lookup_Ac(ri->get_key());

        //NOTE may need to add extra constraint here (for monotonicity?)
        //it's either not the annotated version or x is not optional
        if(!is_xoptional) {
            if(color_v!=NOT_DOMINATED && Aj <= Ai) {
                ri->update_Ac(Aj);
                ri->set_childl_table_ind(rj_tabind);
            }
        } else { //annotated version and x IS an optional vertex.
            if(Aj <= Ai) {
                ri->update_Ac(Aj);
                ri->set_childl_table_ind(rj_tabind);
            }
        }
    }

    return parent_table;
}


Table* update_join_table(Table* rightchildk, Table* leftchildj,
                            Set* optional_verts, Variant variant) {
    /*
     * j is left
     * k is right
     *
     * Initializes parent table.
     *
     * NOTE can combine w/ other join fun
     */
    Table* parent_table = new Table();
    parent_table->vertices = rightchildk->vertices; //NOTE could probably do in the loop

    int tab_size = rightchildk->get_table_size();
    for(int i=0; i<tab_size; i++) {  //4^ni
        Row* row_j = leftchildj->get_row(i);
        Row* row_k = parent_table->create_row(rightchildk->get_row(i), -1);

        int* min_vals = minAi_c(rightchildk, leftchildj, optional_verts, row_k, row_j);
        row_k->update_Ac(min_vals[0]);
        row_k->set_childl_table_ind(min_vals[1]);
        row_k->set_childr_table_ind(min_vals[2]);

        delete[] min_vals;
    }
    return parent_table;
}


//------ Merge Table Versions (reuse a child table, and convert it into the parent)
void merge_introduce_table(Graph* graph, Table* child_table, Set* child_bag, Set* optional_verts,
                            int v, Variant variant) {
    /*
     * Updates the child table to be the parent.
     */
    child_table->vertices.push_back(v);
    Set* neighbors_v = child_bag->set_intersection(graph->neighbors(v));
    int table_size = child_table->get_table_size();

    for(int i=0; i<table_size; i++) {           //3^ni
        Row* r_update = child_table->get_row(i);
        //r_update->set_childl_table_ind(i);

        Row* r2 = child_table->create_row(r_update, DOMINATED);
        Row* r3 = child_table->create_row(r_update, NOT_DOMINATED);

        //r_update: x is IN_DOMSET=1
        intro_vert_indomset_update(graph, child_table, neighbors_v, r_update, v, variant);
        child_table->update_row_add(r_update, IN_DOMSET); //must add after intro vert function call

        //r2: x is DOMINATED=0
        intro_vert_dominated_update(graph, child_table, neighbors_v, optional_verts,
                                    r2, r_update, v, variant);

        //r3: x is NOT_DOMINATED=3  -- do nothing
    }
    delete neighbors_v;
}


void merge_forget_table(Table* child_table, Set* optional_verts, int v, Variant variant) {
    /*
     * Updates the child_table to be the parent.
     */
    int v_index = child_table->get_vertex_col_index(v);      //k
    int table_size = child_table->get_table_size();
    bool is_xoptional = optional_verts->contains(v);  // in optional set or not?

    //delete v from verts vec. k
    child_table->vertices.erase(child_table->vertices.begin()+v_index);

    int curr_index = 0;
    for(int j=0; j<table_size; j++) {                       //3^ni
        Row* row_child = child_table->get_row(curr_index);

        int color_v = row_child->coloring[v_index];
        int Aj = row_child->get_Ac();
        int rj_tabind = child_table->lookup_table_index(row_child->get_key());

        child_table->table_lookups_remove(row_child->get_key());
        row_child->remove_from_coloring(v_index);           //k

        int par_row_ind;
        Row* row_par;
        if(child_table->table_lookups_contains(row_child->get_key())) {
            child_table->delete_row(curr_index);
            curr_index--;

            // gets the actual parent row which has been added previously
            row_par = child_table->lookup_row(row_child->get_key());
        } else {  //row needs readding to table lookups
            par_row_ind = curr_index;
            row_par = child_table->get_row(par_row_ind);

            //add to table lookups.
            child_table->table_lookups_insert(row_par->get_key(), par_row_ind);
        }

        int Ai = child_table->lookup_Ac(row_par->get_key());

        //NOTE may need to add extra constraint here (for monotonicity?)
        //it's either not the annotated version or x is not optional
        if(!is_xoptional) {
            if(color_v!=NOT_DOMINATED && Aj <= Ai) {
                row_par->update_Ac(Aj);
                //row_par->set_childl_table_ind(rj_tabind);
            }
        } else { //annotated version and x IS an optional vertex.
            if(Aj <= Ai) {
                row_par->update_Ac(Aj);
                //row_par->set_childl_table_ind(rj_tabind);
            }
        }

        curr_index++;
    }
}


void merge_join_table(Table* rightchildk, Table* leftchildj,
                       Set* optional_verts, Variant variant) {
    /*
     * j is left
     * k is right  -- calling table
     *
     * i is parent table
     *
     * Updates the rightchild table to be the parent.
     */
    int tab_size = rightchildk->get_table_size();
    for(int i=0; i<tab_size; i++) {  //4^ni
        Row* row_j = leftchildj->get_row(i);
        Row* row_k = rightchildk->get_row(i);

        int* minAi = minAi_c(rightchildk, leftchildj, optional_verts, row_k, row_j);
        row_k->update_Ac(minAi[0]);
        row_k->set_childl_table_ind(minAi[1]);
        row_k->set_childr_table_ind(minAi[2]);

        delete[] minAi;
    }
}



//--- Table reduction helpers
bool is_exclusive_to_singlebag(std::vector<Set*>& bags, Set* reduced_bags,
                            Set* removed_elems, Set* removed_nodes) {
    /* Rule 1: if bag Xi contains at least one element that is not in any other bag,
     * and the element is not contained in the removed elemens set,
     * put Xi in solution and add its elements to the removed elemens set.
     *
     * returns true if something about the solution has been changed by rule.
     * false if input remains exactly the same.
     */
    bool rule1=false;
    for(int i=0; i<bags.size(); i++) {
        bool add_bag=false;
        Set* curr_bag = bags[i];

        //not in soln and hasnt been removed from tree
        if(!reduced_bags->contains(i) && !removed_nodes->contains(i)) {
            for(auto it=curr_bag->begin(); it!=curr_bag->end(); it++) {
                bool unique_v=true;
                int v=*it;

                if(!removed_elems->contains(v)) {
                    for(int j=0; j<bags.size(); j++) {
                        Set* b = bags[j];
                        if(i!=j && !removed_nodes->contains(j) && b->contains(v)) unique_v=false;
                    }

                    if(unique_v) add_bag=true;
                }
            }
        }
        if(add_bag) {
            reduced_bags->insert(i);
            removed_elems->add_all(curr_bag);
            rule1=true;  //added at least on bag
        }
    }
    return rule1;
}


bool is_child_subset(std::vector<Set*>& bags, std::vector<po_bag>& postorder,
                     Set* removed_elems, Set* removed_nodes) {
    /* Rule 2: If there is an edge e = {i, j} ∈ F such that
     *          (1) Xi ⊆ Xj and
     *          (2) w(Xi) ≥ w(Xj)
     * then connect each neighboring node of i (except for j) with j
     * and delete i together with its incident edges from T.
     */
    bool rule2=true;
    bool removed_elem=false;
    for(int i=0; i<postorder.size(); i++) {
        po_bag po = postorder[i];

        int xi_index = po.parent_bag_index;
        int xj_index = po.bag_index;  //current

        //child is xj, par is xi
        //curr po's parent is not the root, and there exists and edge between par and child.
        if(xi_index!=-1 && xi_index!=-2 && !removed_elem) {
            bool subset = is_special_subset(bags[xi_index], bags[xj_index], removed_elems);

            if(subset && pow(3, bags[xi_index]->size()) >= pow(3, bags[xj_index]->size())) {
                //then connect each nbr node of i (except j) w/ j and delete i from T.
                removed_nodes->insert(po.parent_bag_index);
                remove_node_from_postack(postorder, postorder[i]); //remove the parent of i
                removed_elem=true;  //only delete one bag per iteration of while loop
            }
        }
    }
    if(!removed_elem) rule2=false; //did not remove an element

    return rule2;
}


bool is_a_implies_b(Graph* graph, std::vector<Set*> &bags, Set* reduced_bags,
                 Set* removed_elems, Set* removed_nodes) {
    /* Rule 3: If for any bag Xi the fact a ∈ Xi implies b ∈ Xi, then
     * remove b from the bags of the tree decomposition.
     *
     * a vertex b_ is only in a bag when a vertex a_ is in the same bag.
     */
    bool rule3=true;
    bool found_v=false;
    Set* verts = graph->get_vertices();  //n time.

    //n^3 time! NOTE look for better way.
    for(auto it=verts->begin(); it!=verts->end(); it++) {
        int a=*it;
        if(!removed_elems->contains(a)) {
            for(auto itt=verts->begin(); itt!=verts->end(); itt++) {
                int b=*itt;
                if(!removed_elems->contains(b)) {
                    if(a!=b){
                        bool implies=true; //presense of a in bag implies presense of b
                        for(int i=0; i<bags.size(); i++) {
                            Set* curr_bag=bags[i];
                            if(!reduced_bags->contains(i) && !removed_nodes->contains(i)) {
                                if(!curr_bag->contains(a) && curr_bag->contains(b)) implies=false;
                            }
                        }
                        if(implies) {
                            removed_elems->insert(b); //remove b from all bags.
                            found_v=true;
                        }
                    }
                }
            }
        }
    } if(!found_v) rule3=false;
    delete verts;

    return rule3;
}


bool is_empty_bag_intersect(std::vector<Set*> &bags, std::vector<po_bag> &postorder,
                         Set* removed_elems) {
    /* Rule 4: If there is an edge e = {i, j} ∈ F such that Xi ∩ Xj = ∅,
     * then remove this edge from F.
     */
    bool rule4=true;
    bool removed_edge=false;
    for(int i=0; i<postorder.size(); i++) {
        po_bag po = postorder[i];

        int xi_index = po.parent_bag_index;
        int xj_index = po.bag_index;  //current

        //child is xj, par is xi
        //curr po not the root, and there is an edge
        if(xi_index!=-1 && xi_index!=-2) {
            Set* intersect = bags[xj_index]->set_intersection(bags[xi_index]);
            bool empty=true;
            for(auto it=intersect->begin(); it!=intersect->end(); it++) {
                int v=*it;
                if(!removed_elems->contains(v)) empty=false;
            }
            delete intersect;

            if(empty) {
                //removes edge between curr po and parent
                remove_edge_from_postack(postorder, postorder[i]);
                removed_edge=true;
            }
        }
    }
    if(!removed_edge) rule4=false;  //no edges removed
    return rule4;
}


bool is_treelike_subcollection(std::vector<Set*> &bags, std::vector<po_bag> &postorder,
                            Set* removed_elems, Set* removed_nodes) {
    /* Rule 5: If for a bag Xi there is a set of nodes
     * j1, j2, . . . , jl ∈ I, such that
     *      (1) Xi ⊆ UNION k=1 to l of X_jk,
     *      (2) w(Xi) ≥  SUM k=1 to l of w(X_jk),
     *      (3) the neighbors N(i) of i form a tree-like subset collection,
     *
     * then remove i from I and connect its neighbors such
     * that the tree-likeness is not violated
     *
     * xi is parent, xj's are children
     */
    bool rule5=true;
    bool removed_elem2=false;
    for(int i=0; i<postorder.size(); i++) {
        po_bag po = postorder[i];

        //here xi is the parent and current po
        int xi_index = po.bag_index;

        if(po.num_children>1) { //not applicable to nodes w/ only one child
            po_bag child1, child2;
            int po_xj1=-3;  //child 1 po index
            int po_xj2=-3;  //child 2 po index

            //get indices of po_bags of the children in postorder vec
            //if the edge was already deleted, the equality wouldnt pass, so it works
            for(int j=0; j<postorder.size(); j++) {
                if(postorder[j].parent_bag_index==xi_index) {
                    if(po_xj1==-3) {
                        po_xj1=j;
                        child1=postorder[j];
                    } else if(po_xj1!=-3 && po_xj2==-3) {
                        po_xj2=j;
                        child2=postorder[j];
                    }
                }
            }

            int bag_index_xj1 = child1.bag_index;
            int bag_index_xj2 = child2.bag_index;

            Set* un = bags[bag_index_xj1]->set_union(bags[bag_index_xj2]);
            //if only one of the children was subset of parent, it would
            //have been caught by second rule (should have at least).
            bool subset = is_special_subset(un, bags[xi_index], removed_elems);

            int w_xi = pow(3, bags[xi_index]->size());
            int w_xj = pow(3, bags[bag_index_xj1]->size()) + pow(3, bags[bag_index_xj2]->size());
            if(subset && w_xi>=w_xj) {
                removed_nodes->insert(po.bag_index);
                remove_node_from_postack(postorder, postorder[po_xj1]); //remove the parent
                removed_elem2=true;  //only delete one bag per iteration of while loop
            }

            delete un;
        }
    }
    if(!removed_elem2) rule5=false; //did not remove an element

    return rule5;
}


void remove_node_from_postack(std::vector<po_bag> &postorder, po_bag &xj) {
    /* child is xj, par is xi --removing xj's parent (xi)
     *
     * must get xi's parent index
     * set xj/child_po's parent index from xi to xi's parent.
     * if xi has more than one child, set the other childs parent to xi's parent
     *
     * delete xi from T
     */
    int po_xi_ind = -1;  //index in po stack
    int xi_par_ind=-1;   // parent of child's parent (child's grandparent)
    int num_xi_child=-1;

    //get the parents to be removed info
    for(int j=0; j<postorder.size(); j++) {
        if(postorder[j].bag_index==xj.parent_bag_index) po_xi_ind=j;
    }
    xi_par_ind = postorder[po_xi_ind].parent_bag_index;
    num_xi_child = postorder[po_xi_ind].num_children;

    postorder.erase(postorder.begin()+po_xi_ind);

    if(num_xi_child>1) {
        xj.num_children++;
        for(int j=0; j<postorder.size(); j++) {
            if(postorder[j].parent_bag_index==xj.parent_bag_index
                && postorder[j].bag_index!=xj.bag_index) {
                postorder[j].parent_bag_index = xj.bag_index;
            }
        }
    }

    //reset parent of xj (curr po) and if xi had more children, their parent too.
    xj.parent_bag_index = xi_par_ind;

}


void remove_edge_from_postack(std::vector<po_bag> &postorder, po_bag &child_po) {
    //removes edge between curr po and parent

    //sets parent of childs num children--
    int po_xi_ind = -1;  //index in po stack
    //get the parent's po index
    for(int j=0; j<postorder.size(); j++) {
        if(postorder[j].bag_index==child_po.parent_bag_index) po_xi_ind=j;
    }
    postorder[po_xi_ind].num_children--;

    //sets child par index to -2 (no parent but not root)
    child_po.parent_bag_index=-2;
}


bool is_special_subset(Set* A, Set* B, Set* excluded) {
    //Is set A\excluded a subset of set B\excluded?
    bool is_sub = false;
    for(auto it=A->begin(); it!=A->end(); it++) {
        int x = *it;

        if(!excluded->contains(x)) {
            if (!B->contains(x)) is_sub = false;
            else is_sub=true;
        }
    }
    return is_sub;
}


//----- Helper functions

int locally_valid_coloring(Graph* graph, Set* optional_verts, Row* row,
                           std::vector<int> &vertices, Variant variant) {
    /*
     * NOTE Could replace vertices vector w/ bag Set
     * A_c is the size of the dominating set of the specific coloring in the row.
     * If its an invalid coloring, return -99;
     *
     * locally invalid coloring:

        there's some vertex s in the specific coloring which is set to 0
        (ie. already dominated),

        AND there is NOT a vertex, which is a neighbor of the vertex s
        (in the same bag), w/ assigned color 1 (that's IN the dom. set)

        - so the coloring is NOT valid. As in its saying that theres a
        vertex which is dominated, but none of its neighbors are in the
        dom set, so it's not true.

        k^2
     */
    //If finding independent dominating set, create new set.
    Set* independent_check;
    if(variant == Variant::Indep_Dom_Set) independent_check = new Set();

    int A_c = 0;
    for(int i=0; i<vertices.size(); i++) {  //at most k
        int xt = vertices[i];
        int coloring_i = row->coloring[i];

        if(coloring_i == IN_DOMSET) {
            A_c++;
            if(variant == Variant::Indep_Dom_Set) independent_check->insert(xt);
        }

        if(coloring_i == DOMINATED) {  //if a vertex is set to DOMINATED
            bool valid = false;
            bool is_xtoptional = optional_verts->contains(xt);  // in optional set or not?
            int num_dominators = 0;

            //xt is not optional, or it's the annotated version and it is optional
            if(!is_xtoptional) {
                num_dominators = get_num_dominators(graph, row, vertices, xt);
                if(num_dominators>=1) valid=true;
            } else valid = true;

            /* perfect domset variant-each dominated vertex
             * can only be dominated by exactly one vert in domset
             * optional or not
             */
            if(variant == Variant::Perf_Dom_Set) {
                // if x is optional, we still need to check if it
                // has more than one dominator.
                if(is_xtoptional) {
                    num_dominators = get_num_dominators(graph, row, vertices, xt);
                }
                if(num_dominators>1) valid=false;
            }
            //------

            if(!valid) return INF;
        }
    }

    if(variant == Variant::Indep_Dom_Set) {
        bool indep = true;
        if(independent_check->size() >= 2) {
            indep = check_independent(graph, independent_check);
        }

        delete independent_check;
        if(!indep) return INF;  // vertices are not independet
    }
    return A_c;
}


int phi(Row* row, Set* neighbors, std::vector<int> &vertices, int introduced_v) {
   /*
    *  φ : {0, ˆ0, 1}^nj → {0, ˆ0, 1}^nj
    * on the set of colorings of Xj.
    * For c =(c1,... ,c_nj ) ∈ {0, ˆ0, 1}^nj , let φ(c):= (c'_1,... ,c'_nj )
    * such that
    *
    * c'_t = 0ˆ if t ∈ {p1,... ,ps} and ct = 0  OR
    * c'_t = c_t otherwise.
    *
    * returns key of the new coloring
    */
    int new_col;
    for(int i=0; i<vertices.size(); i++) {
        int curr_v = vertices[i];
        if(curr_v != introduced_v) {
            int curr_v_c = row->coloring[i];

            // the introduce vertex has a neighbor w. color set to DOMINATED.
            if(curr_v_c==DOMINATED && neighbors->contains(curr_v)) {
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


int* minAi_c(Table* childk, Table* childj, Set* optional_verts, Row* row_k, Row* row_j) {
    /*
     * NOTE: this is where storing the vertices vector is important. Look into not storing
     * it in table.
     *
     */

    //First, find all possible c' and c'' which divide c.
    std::vector<int> c_prime;
    std::vector<int> c_primeprime;

    int num_ones = 0;   //number of 1's (IN_DOMSET's) in the coloring.
    for(int i=0; i<row_k->coloring.size(); i++) { //k
        int c_t = row_k->coloring[i];

        int xt = childk->vertices[i];
        bool is_xtoptional = optional_verts->contains(xt);  // xt in optional set or not?

        if(c_t == IN_DOMSET) num_ones++;

        int cprime_size = c_prime.size();
        int cprimeprime_size = c_primeprime.size();

        if(cprime_size != cprimeprime_size) printf("ERROR: in divide().\n");

        //is_xoptional will only be true is we want the annotated version.
        if(c_t==IN_DOMSET || c_t==NOT_DOMINATED || is_xtoptional) {
            if(cprime_size==0 && cprimeprime_size==0) {
                int cprime0 = c_t;
                c_prime.push_back(cprime0);    //c_t' == c_t

                int cprimeprime0 = c_t;
                c_primeprime.push_back(cprimeprime0);   //c_t'' == c_t
            } else {
                for(int t=0; t<cprime_size; t++) {
                    //c_t' == c_t
                    c_prime[t] = c_prime[t]*10+c_t; //append c_t to end of int 'string'

                    //c_t'' == c_t
                    c_primeprime[t] = c_primeprime[t]*10+c_t;
                }
            }
        } else if (c_t == DOMINATED ) {
            if(cprime_size==0 && cprimeprime_size==0) {
                //---- c'_t = 0
                int cprime0 = DOMINATED;
                c_prime.push_back(cprime0);

                //c''_t = 0hat
                int cprimeprime0 = NOT_DOMINATED;
                c_primeprime.push_back(cprimeprime0);

                //---- c'_t = 0hat
                int cprime1 = NOT_DOMINATED;
                c_prime.push_back(cprime1);

                //c''_t = 0
                int cprimeprime1 = DOMINATED;
                c_primeprime.push_back(cprimeprime1);
            } else {
                //NOTE combine these for loops
                for(int t=0; t<cprime_size; t++) {
                    //1. make copy
                    int cprime_t = c_prime[t];

                    //2. update current
                    c_prime[t] = cprime_t*10+DOMINATED;

                    //3. add to second pair
                    cprime_t = cprime_t*10+NOT_DOMINATED;

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

    //Now, find Ai(c) ← min{Aj(c')+ Ak(c'') − #1(c) | c' and c'' divide c }
    int min_Ai = INF;
    int A_jcprime_ind_final = -1;
    int A_kprimeprime_ind_final = -1;

    for(int i=0; i<c_prime.size(); i++) {
        int cprime_key = c_prime[i];
        int cprimeprime_key = c_primeprime[i];

        int A_jcprime_ind = childj->lookup_table_index(cprime_key);
        int A_jcprime = childj->lookup_Ac(cprime_key);

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

    // j is left
    // k is right
    int *min_vals = new int[3];
    min_vals[0] = min_Ai;
    min_vals[1] = A_jcprime_ind_final;
    min_vals[2] = A_kprimeprime_ind_final;

    return min_vals;
}


void intro_vert_indomset_update(Graph* graph, Table* child_table,
                                Set* neighbors_v, Row* row,
                                int v, Variant variant) {
    /* Introducing a vertex w/ coloring set to IN_DOMSET.
     */
    int new_col_key = phi(row, neighbors_v, child_table->vertices, v);   //k
    int A_phi;
    bool is_perf=true;
    //Independent variant
    if(variant == Variant::Indep_Dom_Set) {
        bool indep = intro_indep_check(graph, child_table->vertices, row->coloring, v);
        if(!indep) A_phi = INF;
        else {
            A_phi = child_table->lookup_Ac(new_col_key);
            if(A_phi!=INF) A_phi++;
        }
    } else if(variant == Variant::Perf_Dom_Set) {
        //gets neighbors of the newly added vertex
        //if a neighbor is not in domset and has more than one neighbor in domset, Ac=inf
        for(auto it=neighbors_v->begin(); it!=neighbors_v->end(); it++) {
            int nb = *it;
            int c_index = child_table->get_vertex_col_index(nb);
            int cp = row->coloring[c_index];
            if(cp==DOMINATED) {
                int num_doms_rnb = get_num_dominators(graph, row, child_table->vertices, nb);
                if(graph->adjacent(nb, v)) num_doms_rnb++;  //v not in the child table vertices vec.
                if(num_doms_rnb > 1) is_perf = false;
            }
        }
        if(!is_perf) A_phi = INF;
        else {
            A_phi = child_table->lookup_Ac(new_col_key);
            if(A_phi!=INF) A_phi++;
        }
    } else {
        A_phi = child_table->lookup_Ac(new_col_key);  //Dom Set variant
        if(A_phi!=INF) A_phi++;
    }
    row->update_Ac(A_phi);
}


void intro_vert_dominated_update(Graph* graph, Table* child_table,
                                 Set* neighbors_v, Set* optional_verts,
                                 Row* r2, Row* r_update, int v, Variant variant) {
    /* Update to the introduce table when the introduced vertex coloring
     * is set to DOMINATED.
     */
    bool is_xoptional = optional_verts->contains(v);  // in optional set or not?
    int num_dominators_r2 = 0;
    bool ru_isperfect = true;

    //it's either not the annotated version or x is not optional
    if(!is_xoptional) {
        //need to check that the introduced vertex w. the coloring of DOMINATED is justified.
        for(auto it=neighbors_v->begin(); it!=neighbors_v->end(); it++) { //*k^2!!
            int index = child_table->get_vertex_col_index(*it);

            //x is set to DOMINATED and has a neighbor IN_DOMSET
            if(r2->coloring[index]==IN_DOMSET) {
                num_dominators_r2++;
            }

            if(variant==Variant::Perf_Dom_Set) {
                if(r_update->coloring[index]==DOMINATED) {
                    //checks if the current neighbor of x (xt) has
                    //more than one dominator now.
                    int xt = child_table->vertices[index];
                    int num_doms_rup_nb = get_num_dominators(graph, r_update,
                                                             child_table->vertices, xt);
                    if(num_doms_rup_nb > 1) ru_isperfect=false;
                }
            }
        }
        if(num_dominators_r2<1) {  //not valid coloring
            r2->update_Ac(INF);
        }
    } else {
        //annotated version and x IS an optional vertex, dont need to justify.
        //Ai(c x {DOMINATED}) <- Aj(c)  i.e. do nothing

        //same checks for Perfect variant though
        if(variant==Variant::Perf_Dom_Set) {
            for(auto it=neighbors_v->begin(); it!=neighbors_v->end(); it++) { //*k^2!!
                int index = child_table->get_vertex_col_index(*it);
                if(r2->coloring[index]==IN_DOMSET)  num_dominators_r2++;
                if(r_update->coloring[index]==DOMINATED) {
                    int xt = child_table->vertices[index];
                    int num_doms_rup_nb = get_num_dominators(graph, r_update,
                                                            child_table->vertices, xt);
                    if(num_doms_rup_nb > 1) ru_isperfect=false;
                }
            }
        }
    }

    if(variant==Variant::Perf_Dom_Set) {
        //More than one dominator for the perfect ds variant not allowed.
        if(num_dominators_r2>1) r2->update_Ac(INF);
        if(!ru_isperfect) r_update->update_Ac(INF);
    }
}


int get_num_dominators(Graph* graph, Row* row, std::vector<int> &vertices, int target_vert) {
    //returns number of vertices dominating target_vert.
    int num_dominators=0;
    for(int j=0; j<vertices.size(); j++) {  //at most k
        int coloring_j = row->coloring[j];

        if(vertices[j]!=target_vert && coloring_j== IN_DOMSET) {
            if(graph->adjacent(vertices[j], target_vert)) {
                num_dominators++;
            }
        }
    }
    return num_dominators;
}


bool intro_indep_check(Graph* graph, std::vector<int> &vertices, std::vector<int> &coloring,
                       int v) {
    Set* independent_check = new Set();
    for(int i=0; i<coloring.size(); i++) {
        int c=coloring[i];
        int xt=vertices[i];

        if(c==IN_DOMSET) {
            independent_check->insert(xt);
        }
    }
    independent_check->insert(v);  //check against the introduced vert.

    bool check = check_independent(graph, independent_check); //k^2
    delete independent_check;

    return check;
}


bool check_independent(Graph* graph, Set* independent_check) {
    /*
     * returns true if none of the vertices in the set are adjacent (are independent),
     * returns false is two vertices are adjacent.
     */
    for(auto it=independent_check->begin(); it!=independent_check->end(); it++) {
        int u = *it;
        for(auto itt=independent_check->begin(); itt!=independent_check->end(); itt++) {
            int v = *itt;
            if(u!=v && graph->adjacent(u, v)) return false; // are adjacent
        }
    }
    return true;
}


//---- For testing
void print_row(Row* row) {
    printf("|");
    for(int k=0; k<row->coloring.size(); k++) {
        printf("  %d  |", row->coloring[k]);
    }
    printf("%15d \t | ", row->get_Ac());
    printf("%15d \t", row->get_childl_table_ind());
    printf("%5d \t", row->get_childr_table_ind());
    printf("\n");
}

void print_table(Table* tab, std::string table_type) {
    printf("\n=====================================================================\n");
    printf("IN_DOMSET=1=1,   DOMINATED=2=0,   NOT_DOMINATED=3=0hat\n");
    printf("            Type: %s  \n", table_type.c_str());
    printf("  vertices \t\t |  A_ci \t |  childl ind  |  childr_ind\n");
    printf("|");
    for(int j=0; j<tab->vertices.size(); j++) {
        printf("  %d  |", tab->vertices[j]);
    }
    printf("\n____________________________________________________________________\n");
    printf("\n");

    for(int j=0; j<tab->get_table_size(); j++) {
        Row* r = tab->get_row(j);
        print_row(r);
    }

    printf("\n=====================================================================\n\n\n");
}

void print_tables(std::vector<Table*> &tables) {
    for(int i=0; i<tables.size(); i++) {
        Table* tab = tables[i];

        std::string type = "All";
        print_table(tab, type);
    }
}

void print_postorder(std::vector<po_bag> postorder) {
    printf("\n------POSTorder\n");

    for(int j=0; j<postorder.size(); j++) {
        //printf("postorder: ind=%d, num_childs=%d\n", po[j].bag_index, po[j].num_children);
        print_pobag(postorder[j]);
    }
    printf("---------\n\n");
}

void print_pobag(po_bag po) {
    printf("bag index=%d\n", po.bag_index);
    printf("num children=%d\n", po.num_children);
    printf("parent bag index=%d\n", po.parent_bag_index);
    printf("\n");
}
