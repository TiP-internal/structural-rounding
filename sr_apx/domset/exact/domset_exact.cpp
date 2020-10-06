
#include "domset_exact.hpp"

#include <cstdio>
#include <cassert>


//--- Calculating Solution
int get_soln_row_index(Table* table, Set* optional_verts){
    //Returns the index of the row in the last table w/ smallest dominating set size.
    int soln_index=-1;
    int min_Ai = INF;

    //finds the row w. smallest solution size.
    for(int i=0; i<table->get_table_size(); i++) {
        Row* row = table->get_row(i);
        
        if(row->get_Ac() < min_Ai) {
            int num_notdom=0;
            for(int j=0; j<table->vertices.size(); j++) {
                int c = row->coloring[j];
                int xt = table->vertices[j];
                //a final coloring cant have any vertices colored as not dominated
                if(c==NOT_DOMINATED && !optional_verts->contains(xt)) num_notdom++;
            }
            
            if(num_notdom==0) {
                min_Ai=row->get_Ac();
                soln_index=i;
            }
        }
    }
    return soln_index;
}


void add_to_solution(Set* dom_set, Row* row, std::vector<int> &vertices) {
    //Adds the vertices which have their coloring set to IN_DOMSET in the given row.
    for(int i=0; i<vertices.size(); i++) {
        int xt = vertices[i];
        int c = row->coloring[i];
        
        if(c==IN_DOMSET) {
            dom_set->insert(xt);
        } 
    }
}


int get_solution(std::vector<Table*> &tables, Set* optional_verts) {
    /*
     * Finds the size of the final solution for the optimization version of DOMSet
     */
    int soln_size=0;
    int tabsize = tables.size();  //shouldnt be more than 2
    Table* tab;
    
    for(int i=0; i<tabsize; i++) {
        tab = tables.back();
        tables.pop_back();
                
        int soln_index = get_soln_row_index(tab, optional_verts);  //prints the tabel lookups
        if(soln_index==-1) printf("ERROR: soln row not found.\n");
        
        Row* final_row = tab->get_row(soln_index);
        soln_size += final_row->get_Ac();

        delete tab;
    }
    return soln_size;
}


int get_solution(std::vector<Table*> &tables, Set* dom_set, Set* optional_verts) {  
    /*
     * This is the second pass of the alg which constructs the soln set.
     * It iterates over tables from top to bottom.
     * 
     * returns number of domset vertices added for this one compoenent
     */
    Table* parent_table = tables.back();
    tables.pop_back();
    
    po_bag parent_pobag = parent_table->get_pobag();
    int soln_index = get_soln_row_index(parent_table, optional_verts);
    
    if(soln_index==-1) printf("ERROR: soln row not found.\n");
    Row* soln_row = parent_table->get_row(soln_index);

    int soln_size = soln_row->get_Ac();
    if(soln_size>=INF) printf("ERROR: in get_solution for constructive version.\n");
    
    add_to_solution(dom_set, soln_row, parent_table->vertices);

    int childl_row_ind = soln_row->get_childl_table_ind();
    int childr_row_ind = soln_row->get_childr_table_ind();

    std::deque<po_bag> unfinished_joins;
    std::deque<int> left_inds;
    int tabsize = tables.size();  //shouldnt be more than 2
    
    for(int i=0; i<tabsize; i++) {
        // goal is to find the correct current row from the child table during each iteration.
        Row* curr_row;  
                
        if(parent_pobag.num_children==2) {        //root table is a join table
            unfinished_joins.push_front(parent_pobag);
            left_inds.push_front(childl_row_ind);
        }
        
        Table* child_table = tables.back();
        tables.pop_back();
        po_bag child_pobag = child_table->get_pobag();  //par_bag_index, bag_index, num_children
       
        if(parent_pobag.bag_index != child_pobag.parent_bag_index) {
            // if the child table's parent is not the current parent
            // then it must must be the left child of an unfinished join table
            if(unfinished_joins[0].bag_index == child_pobag.parent_bag_index) {
                parent_pobag = unfinished_joins[0];  //get most recently added
                unfinished_joins.pop_front();
                
                childl_row_ind = left_inds[0];
                left_inds.pop_front();
                
                curr_row = child_table->get_row(childl_row_ind);
            } else {
                printf("ERROR in getting solution\n");
            }
        } else {
            // now that the current parent has been evaluated, 
            // reassign child table as the parent
            curr_row = child_table->get_row(childr_row_ind); //there will be a left ind if its not a leaf
        }
        add_to_solution(dom_set, curr_row, child_table->vertices);
        
        delete parent_table;
        parent_table = child_table; 
        parent_pobag = child_pobag;
            
        childl_row_ind = curr_row->get_childl_table_ind();
        childr_row_ind = curr_row->get_childr_table_ind();
    }
    
    if(soln_index==-1) printf("ERROR: soln row not found.\n");
    return soln_size;
}


int calculate(Graph* graph, TreeDecomp* decomp, Set* dom_set,
              Set* optional_verts, Variant variant, bool construct_soln) {
    /* Calculates the minimum dominating set given a tree decomp.
     * 
     * if construct_soln is true, the solution set is constructed, 
     * else its not
     *
     * Must calculate for each component in the decomp.
     * Calculates annotated or regular version of the dom set.
     */
    std::vector<Table*> tables; 
    std::vector<po_bag> postorder = decomp->get_post_order();
    
    // Set* anchors = treedecomp_reduction(graph, decomp->components_bags, postorder); 
    Set* anchors = new Set();
    if(construct_soln) for(int i=0; i<decomp->components_bags.size(); i++) anchors->insert(i);
    
    int solutionsize=0;
    for(int i=0; i<postorder.size(); i++) { 
        int bag_index = postorder[i].bag_index;
        int par_bag_index = postorder[i].parent_bag_index;
        
        if(decomp->components_bags[bag_index]->size() != 0) {
            calculate_tables(graph, decomp->components_bags,
                             postorder, tables, 
                             optional_verts,
                             anchors, variant, i);
        }
        
        if(par_bag_index==-1 || decomp->components_bags[par_bag_index]->size() == 0) {  
            //-----------------------empty root bag or the actual root has been reached
            //calculate solution up until this point, and discard current tables       
            if(tables.size()>0) {
                if(construct_soln) {
                    solutionsize += get_solution(tables, dom_set, optional_verts); 
                }
                else {
                    solutionsize += get_solution(tables, optional_verts);
                }
            }
        } 
    }
    delete anchors;
    return solutionsize;
}


int get_left_join_child_tabind(Set* anchor_tables, 
                               std::vector<Table*> &tables,
                               int current_index) {
    for(int j=tables.size()-2; j>=0; j--) {  //-2 since tables.size()-1 is the right child
        //if(!anchor_tables->contains(ind) && tab_ind_childr==-99) {
        po_bag check = tables[j]->get_pobag();
        if(check.parent_bag_index==current_index) {
            return j;  //left child table index
        }
    }
    return -1;
}


void calculate_tables(Graph* graph, std::vector<Set*> &bags,
                      std::vector<po_bag> &postorder, 
                      std::vector<Table*> &tables,
                      Set* optional_verts, Set* anchor_tables, 
                      Variant variant, int po_index) {
    /*
     * Dynamic programming algorithm, for dominating set on graphs
     * w/ bounded treewidth.
     *
     * Constructive version and does not reuse (all) tables. To save memory, it
     * will eventually only store anchor tables.
     *
     * For NICE tree decompositions. 
     * For both annotated and regular dominating set, 
     * and annotated/regular versions of independent dominating set 
     * and perfect dominating set. 
     */
    po_bag pob_current = postorder[po_index];
    
    Table* table = nullptr;
    if(pob_current.num_children==2) {          //-----------------------JOIN bag
        //---------------Get right child table index in tables[]
        // right child will always be the most recently added
        int tab_ind_childr = tables.size()-1;  
        Table* right_child_table = tables[tab_ind_childr];
        po_bag pob_rightchild = right_child_table->get_pobag();
        if(pob_rightchild.parent_bag_index != pob_current.bag_index) {
            printf("ERROR: parent (right) child dont align.\n");
        }
        int bag_ind_childr = pob_rightchild.bag_index;
        //----------------
        
        //---------------Get left child table index in tables[]
        int tab_ind_childl = get_left_join_child_tabind(anchor_tables, tables, pob_current.bag_index);
        if(tab_ind_childl==-1) {
            printf("ERROR: parent (left) child not in tables vector.\n");
        }
        Table* left_child_table = tables[tab_ind_childl];
        po_bag pob_leftchild = left_child_table->get_pobag();
        if(pob_leftchild.parent_bag_index != pob_current.bag_index) {
            printf("ERROR: parent (left) child dont align.\n");
        }
        int bag_ind_childl = pob_leftchild.bag_index;
        //----------------
        
        //printf("bag_ind_childl=%d, bag_ind_childr=%d\n", bag_ind_childl, bag_ind_childr);
        if(left_child_table->get_table_size() != right_child_table->get_table_size()) {
            printf("ERROR in calculate_tables: join children have different dimension\n");
        }

        //update table here
        //join_table(tab1, tab2, ...) -- tab2 is merged into tab1 to become the parent of tab1 & tab2
        if(anchor_tables->contains(bag_ind_childr) && anchor_tables->contains(bag_ind_childl)) {
            //if both children are anchor tables, create new table and add to tables  NOTE HERE
            bool merge = false;
            table = join_table(right_child_table, left_child_table, optional_verts, pob_current, merge);
        } else {
            //if neither are anchors, merge right child into parent table, delete left
            bool merge = true;
            table = join_table(right_child_table, left_child_table, optional_verts, pob_current, merge);
            
            delete left_child_table;
            Table* left_child_table=nullptr;
            tables.pop_back();
            tables.erase(tables.begin()+tab_ind_childl); // deletes the left child table.
        }
    } else if(pob_current.num_children==1) {     //-------------------either INTRODUCE or FORGET bag
        po_bag pob_child = postorder[po_index-1];
        int child_bag_index = pob_child.bag_index;
        int child_bag_table_index = tables.size()-1;  //child is the most recently added to end of tables vec

        Set* parent_bag = bags[pob_current.bag_index];
        Set* child_bag = bags[child_bag_index];

        Table* child_table = tables[child_bag_table_index];
        
        int v;
        if(parent_bag->size() > child_bag->size()) {        // introduce node
            // gets the introduced vertex
            for(auto it=parent_bag->begin(); it!=parent_bag->end(); it++) {
                int u = *it;
                if(!child_bag->contains(u)) v = u;
            }
            
            bool merge;
            if(anchor_tables->contains(child_bag_index)) {
                //create a new table if the child is an anchor
                merge = false;
            } else {
                //child table is now converted into the current par table
                merge = true;
                tables.pop_back();
            }
            table = intro_table(graph, child_table, child_bag, optional_verts, pob_current, variant, v, merge); 
            
        } else if(parent_bag->size() < child_bag->size()) { // forget node
            // gets the forgotten vertex
            for(auto it=child_bag->begin(); it!=child_bag->end(); it++) {
                int u = *it;
                if(!parent_bag->contains(u)) v = u;
            }

            bool merge;
            if(anchor_tables->contains(child_bag_index)) {
                //create a new table if the child is an anchor
                merge=false;
            } else {
                merge=true;
                tables.pop_back();
            }
            table = forget_table(child_table, optional_verts, pob_current, variant, v, merge);
        } else {
            printf("ERROR in calculate_tables: parent and child should not have same size.\n");
        }
    } else if(pob_current.num_children==0){    //----------------------------LEAF bag
        table = initialize_leaf_table(graph, bags[pob_current.bag_index], optional_verts, pob_current, variant);
    } else {
            printf("ERROR in calculate_tables: wrong number of children in nice decomp.\n");
    }
    
    //------ now Add/and or update the pointer in the tables vector.
    tables.push_back(table); 
}


Table* init_parent_table(Table* child1_table, bool merge) {
    Table* par_table = nullptr;
    if(!merge) { //creating new table  (constructive version)
        par_table = new Table();
        par_table->vertices = child1_table->vertices; // k time
    }else if(merge) { //merging    (decision / or child1 (and/or child2) is not an anchor table)
        par_table = child1_table;  //assign the pointer par_table to have the same value as child1_table
    }
    return par_table;
}

Table* initialize_leaf_table(Graph* graph, Set* bag, Set* optional_verts, po_bag pob_current, Variant variant) {
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
    
    table->set_pobag(pob_current);
    return table;
}

Table* join_table(Table* child1_table, Table* child2_table, Set* optional_verts, po_bag pob_current, bool merge) {
    /*
     * If merge is true, the par_table is set to the child_table and is not a nullptr 
     *      --then par_table has already been set to the child1_table. 
     * 
     * if merge if false, set par_table = new Table();, and parent_table->vertices = child1_table->vertices; 
     */
    Table* par_table = init_parent_table(child1_table, merge);
    
    if(child2_table == nullptr) printf("ERROR: child2_table should NOT be a nullptr\n");
    if(child1_table == nullptr) printf("ERROR: child1_table should NOT be a nullptr\n");
        
    int tab_size = child2_table->get_table_size();
    for(int i=0; i<tab_size; i++) {  //4^ni        
        Row* row_k = nullptr;
        if(!merge) {  //creating new table  (constructive version)
            if(child1_table == nullptr) printf("ERROR: child1_table should NOT be a nullptr at this point\n");
            row_k = par_table->create_row(child1_table->get_row(i), -1); //declares new
        } else {      //merging
            row_k = par_table->get_row(i);  //row_k differs between two versions
        }
        
        int* min_vals = minAi_c(child1_table, child2_table, optional_verts, row_k);
        row_k->update_Ac(min_vals[0]);
        row_k->set_childl_table_ind(min_vals[1]);
        row_k->set_childr_table_ind(min_vals[2]);  //merge into row_k in par_table

        delete[] min_vals;
    }
    
    if(merge) child1_table=nullptr;
    par_table->set_pobag(pob_current);
    return par_table;
}

Table* intro_table(Graph* graph, Table* child_table, Set* child_bag, Set* optional_verts,
                   po_bag pob_current, Variant variant, int v, bool merge) {
    /* parent and child same when merging
     * Combined function for both merging and initializing
     */
    Table* par_table = init_parent_table(child_table, merge);
    
    par_table->vertices.push_back(v);
    Set* neighbors_v = child_bag->set_intersection(graph->neighbors(v));
    int table_size = child_table->get_table_size();
    
    for(int i=0; i<table_size; i++) {
        Row* r_update;  //rows in the child_table that we are duplicating and copying
        if(!merge) {
            r_update = new Row(child_table->get_row(i));  //creates copy of row
            par_table->insert_row(r_update);
        } else {
            r_update = par_table->get_row(i);
        }
        r_update->set_childl_table_ind(-1); 
        r_update->set_childr_table_ind(i);
        
        Row* r2 = par_table->create_row(r_update, DOMINATED);
        Row* r3 = par_table->create_row(r_update, NOT_DOMINATED);
        
        //r_update: x is IN_DOMSET=1
        intro_vert_indomset_update(graph, child_table, neighbors_v, r_update, v, variant);
        par_table->update_row_add(r_update, IN_DOMSET); //must add after intro vert function call

        //r2: x is DOMINATED=0
        intro_vert_dominated_update(graph, child_table, neighbors_v, optional_verts,
                                    r2, r_update, v, variant);

        //r3: x is NOT_DOMINATED=3  -- do nothing
    }
    
    if(merge) child_table=nullptr;
    par_table->set_pobag(pob_current);
    delete neighbors_v;
    return par_table;
}

Table* forget_table(Table* child_table, Set* optional_verts, po_bag pob_current, 
                    Variant variant, int v, bool merge) {
    /*
     * if merge is false: Initializes parent table. 
     * if merge is true: Updates the child_table to be the parent.
     */    
    Table* par_table = init_parent_table(child_table, merge);
    int v_index = par_table->get_vertex_col_index(v);  //k
    
    int table_size = child_table->get_table_size();
    bool is_xoptional = optional_verts->contains(v);  // v in optional set or not?
    
    par_table->vertices.erase(par_table->vertices.begin()+v_index); //delete v from verts vec. k
    
    int curr_index = 0;
    for(int j=0; j<table_size; j++) {                       //3^ni        
        Row* rowj=nullptr; 
        if(!merge) {
            rowj = new Row(child_table->get_row(j));  //copy construct k
        } else {
            rowj = child_table->get_row(curr_index);
        }
        
        if(rowj==nullptr) printf("null ptr\n");

        int color_v = rowj->coloring[v_index];
        int Aj = rowj->get_Ac();
        int rj_tabind = child_table->lookup_table_index(rowj->get_key());
        
        if(merge) par_table->table_lookups_remove(rowj->get_key());
        rowj->remove_from_coloring(v_index);  //k
        
        Row* rowi;
        if(!merge) {
            //if(par_table->lookup_table_index(rowj->get_key()) == -1) { //not inserted yet
            if(!par_table->table_lookups_contains(rowj->get_key())) {
                par_table->insert_row(rowj);
            }
            rowi = par_table->lookup_row(rowj->get_key());
        } else {   
            int par_row_ind;
            if(par_table->table_lookups_contains(rowj->get_key())) {  //already inserted in the child
                par_table->delete_row(curr_index);
                curr_index--;

                // gets the actual parent row which has been added previously
                rowi = par_table->lookup_row(rowj->get_key());
            } else {  //row needs readding to table lookups
                par_row_ind = curr_index;
                rowi = par_table->get_row(par_row_ind);

                //add to table lookups.
                par_table->table_lookups_insert(rowi->get_key(), par_row_ind);
            }
        }
        int Ai = par_table->lookup_Ac(rowi->get_key());
        
        //NOTE may need to add extra constraint here (for monotonicity?)
        //it's either not the annotated version or x is not optional
        if(!is_xoptional) {
            if(color_v!=NOT_DOMINATED && Aj <= Ai) {
                rowi->update_Ac(Aj);
                rowi->set_childr_table_ind(rj_tabind);
                rowi->set_childl_table_ind(-1);
            }
        } else { //annotated version and x IS an optional vertex.
            if(Aj <= Ai) {
                rowi->update_Ac(Aj);
                rowi->set_childr_table_ind(rj_tabind);
                rowi->set_childl_table_ind(-1);
            }
        }
        curr_index++;
    }
    
    if(merge) child_table=nullptr;
    par_table->set_pobag(pob_current);
    return par_table;
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


int* minAi_c(Table* childk, Table* childj, Set* optional_verts, Row* row_k) {
    /*
     * This is where storing the vertices vector is important. Look into not storing
     * it in table.
     *
     */
    //First, find all possible c' and c'' which divide c.
    std::vector<std::vector<int>> c_prime;
    std::vector<int> c_primeprime;  //colorings for child k table
    
    bool flip=false;
    for(int i=0; i<childj->vertices.size(); i++) {
        int xt = childj->vertices[i];
        int other = childk->vertices[i];
        if(xt != other) flip=true;
    }
    
    int num_ones = 0;   //number of 1's (IN_DOMSET's) in the coloring.
    for(int i=0; i<row_k->coloring.size(); i++) { //k
        int c_t = row_k->coloring[i];

        int xt = childk->vertices[i];
        bool is_xtoptional = optional_verts->contains(xt);  // xt in optional set or not?

        if(c_t == IN_DOMSET) num_ones++;

        int cprime_size = c_prime.size();
        int cprimeprime_size = c_primeprime.size();

        if(cprime_size != cprimeprime_size) printf("ERROR: in func minAi_c() -- divide.\n");

        //is_xoptional will only be true is we want the annotated version.
        if(c_t==IN_DOMSET || c_t==NOT_DOMINATED || is_xtoptional) {
            if(cprime_size==0 && cprimeprime_size==0) {
                //int cprime0 = c_t;
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
                //---- c'_t = 0
                std::vector<int> cprime0;
                cprime0.push_back(DOMINATED);
                c_prime.push_back(cprime0);

                //c''_t = 0hat
                int cprimeprime0 = NOT_DOMINATED;
                c_primeprime.push_back(cprimeprime0);

                //---- c'_t = 0hat
                std::vector<int> cprime1;
                cprime1.push_back(NOT_DOMINATED);
                c_prime.push_back(cprime1);

                //c''_t = 0
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
    //Now, find Ai(c) ← min{Aj(c')+ Ak(c'') − #1(c) | c' and c'' divide c }
    int min_Ai = INF;
    int A_jcprime_ind_final = -1;
    int A_kprimeprime_ind_final = -1;

    for(int i=0; i<c_prime.size(); i++) {
        //int cprime_key = c_prime[i];
        std::vector<int> cprime_key = c_prime[i];
        int cprimeprime_key = c_primeprime[i];

        int cprime_key_val=0;
        //the keys now need to be flipped so they align with js coloring order
        if(flip) {
            int flipped_cprime_key = flip_coloring(childj, childk->vertices, cprime_key);
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

    // j is left
    // k is right
    int *min_vals = new int[3];
    min_vals[0] = min_Ai;
    min_vals[1] = A_jcprime_ind_final;
    min_vals[2] = A_kprimeprime_ind_final;

    return min_vals;
}


int flip_coloring(Table* childj, std::vector<int> &current_vertorder, std::vector<int> &cprime_key) {
    int flipped_key=-1;
    
    for(int i=0; i<childj->vertices.size(); i++) {  //left
        for(int k=0; k<current_vertorder.size(); k++) { //right
            int xt = childj->vertices[i];
            int other = current_vertorder[k];
            if(xt == other && i==0) flipped_key=cprime_key[k];
            else if(xt==other) flipped_key = flipped_key*10+cprime_key[k];
        }
    }
    
    return flipped_key;
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
        if(A_phi!=INF) A_phi++;  // increase by one
    }
    
    int index = child_table->lookup_table_index(new_col_key);
    row->set_childl_table_ind(-1); 
    row->set_childr_table_ind(index);
    row->update_Ac(A_phi);
}


void intro_vert_dominated_update(Graph* graph, Table* child_table,
                                 Set* neighbors_v, Set* optional_verts,
                                 Row* r2, Row* r_update, int v, Variant variant) {
    /* Update to the introduce table when the introduced vertex coloring
     * is set to DOMINATED.
     * 
     * NOTE get rid of r_update in this function
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
        if(!ru_isperfect) r2->update_Ac(INF);  //NOTE r2 or r_update?
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

