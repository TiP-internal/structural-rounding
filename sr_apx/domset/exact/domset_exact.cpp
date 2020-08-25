
#include "domset_exact.hpp"

#include <cstdio>


//---for testing
bool is_domset(Graph* graph, Set* domset) {    
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int v=*it;
        bool adjacent = false;
        
        for(auto ids=domset->begin(); ids!=domset->end(); ids++) {
            int u=*ids;
            if(graph->adjacent(v, u)) {
                adjacent=true;
            }
        }
        if(!adjacent && !domset->contains(v)) return false;
    }
    return true;
}

bool is_ann_domset(Graph* graph, Set* domset, Set* NX) {
    /*
     * For the annotated version of DS.  
     * 
     * B=V\N[X] MUST be dominated.
     * 
     * N[X] CAN be dominated, but its optional. 
     */
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int v=*it;
        bool adjacent = false;
        for(auto ids=domset->begin(); ids!=domset->end(); ids++) {
            int u=*ids;
            if(graph->adjacent(v, u)) {
                adjacent=true;
            } else if (NX->contains(v)) {  //Doenst need to be dominated.
                adjacent=true;
            }
        }
        if(!adjacent && !domset->contains(v)) return false;
    }
    return true;
}

bool is_indepen_domset(Graph* graph, Set* domset) {
    //is domset    
    bool is_ds = is_domset(graph, domset);
    if(!is_ds) return false;
    
    //is independent
    for(auto ids=domset->begin(); ids!=domset->end(); ids++) {
        int u=*ids;
        for(auto idss=domset->begin(); idss!=domset->end(); idss++) {
            int v=*idss;
            if(graph->adjacent(v, u)) {
                return false;
            }
        }
    }
    return true;
}

bool is_indepen_ann_domset(Graph* graph, Set* domset, Set* NX) {
    //is ann domset
    bool is_ads = is_ann_domset(graph, domset, NX);
    if(!is_ads) return false;
    
    //is independent
    for(auto ids=domset->begin(); ids!=domset->end(); ids++) {
        int u=*ids;
        for(auto idss=domset->begin(); idss!=domset->end(); idss++) {
            int v=*idss;
            if(graph->adjacent(v, u)) {
                return false;
            }
        }
    }
    return true;
}

bool is_perf_domset(Graph* graph, Set* domset) {
    //is domset    
    bool is_ds = is_domset(graph, domset);
    if(!is_ds) return false;
    
    //is perfect
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int u=*it;
        
        if(!domset->contains(u)) {
            int count = 0;
            for(auto itt=domset->begin(); itt!=domset->end(); itt++) {
                int v=*it;
                if(u!=v && domset->contains(v)) count++;
                if(count > 1) return false;
            }
        }
    }
    return true;
}

bool is_per_ann_domset(Graph* graph, Set* domset, Set* NX) {
    //is ann domset    
    bool is_ads = is_ann_domset(graph, domset, NX);
    if(!is_ads) return false;
    
    //is perfect
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int u=*it;
        
        if(!domset->contains(u)) {
            int count = 0;
            for(auto itt=domset->begin(); itt!=domset->end(); itt++) {
                int v=*it;
                if(u!=v && domset->contains(v)) count++;
                if(count > 1) return false;
            }
        }
    }
    return true;
}


void print_row(Row* row) {
    printf("|");
    for(int k=0; k<row->coloring.size(); k++) {
        printf("  %d  |", row->coloring[k]);
    }
    printf("%15d \t | ", row->get_Ac());
    for(auto it=row->domset_verts->begin(); it!=row->domset_verts->end(); it++) {
        printf(" %d,", *it);
    }
    printf("\n");
}

void print_table(Table* tab, std::string table_type) {
    printf("\n===========================================\n");
    printf("IN_DOMSET=1=1,   DOMINATED=2=0,   NOT_DOMINATED=3=0hat\n");
    printf("            Type: %s  \n", table_type.c_str());
    printf("  vertices \t\t |  A_ci \t | Soln. Set \n");
    printf("|");
    for(int j=0; j<tab->vertices.size(); j++) {
        printf("  %d  |", tab->vertices[j]);
    }
    printf("\n___________________________________________\n");
    printf("\n");
    
    for(int j=0; j<tab->get_table_size(); j++) {
        Row* r = tab->get_row(j);
        print_row(r);
    }
    
    printf("\n===========================================\n\n\n");
}

void print_tables(std::vector<Table*> tables) {
    for(int i=0; i<tables.size(); i++) {
        Table* tab = tables[i];
        
        std::string type = "All";
        print_table(tab, type);
    }
}


//---

Set* get_solution(Table* table) {
    /*
     * Constructs the final solution. 
     * Adds the vertices in the final table coloring to the soln set.
     * These are vertices which have not been 'forgotten' vertices yet. 
     */
    printf("Getting solution\n");
    
    int soln_index=-1;
    int min_Ai = INF;
    
    //finds the row w. smallest solution size.
    for(int i=0; i<table->get_table_size(); i++) {
        Row* row = table->get_row(i);
        printf("row->Ai=%d\n", row->get_Ac());
        
        if(row->get_Ac() < min_Ai) {
            min_Ai=row->get_Ac();
            soln_index=i;
        }
    }
    
    if(soln_index==-1) printf("ERROR: soln row not found.\n");
        
    //Add verts (set as IN_DOMSET) in the coloring to the soln set. 
    Row* soln_row = table->get_row(soln_index);
    for(int i=0; i<soln_row->coloring.size(); i++) {
        if(soln_row->coloring[i]==IN_DOMSET) {
            soln_row->domset_verts->insert(table->vertices[i]);         //NOTE for testing
        }
    }
    
    return soln_row->domset_verts;
}


std::vector<Set*> calc_domset(Graph* graph, TreeDecomp* decomp, 
                              Set* optional_verts, Variant variant) {
    /*
     * Calculates the minimum dominating set given a tree decomp.
     * Must calculate for each component in the decomp. 
     * 
     * Calculates annotated and regular version (separately) of the Dominating
     * Set problem. 
    */
    std::vector<std::vector<po_bag>> postorder = decomp->get_post_order();

    Set* dom_set;
    std::vector<Set*> component_domsets;
    for(int j=0; j<decomp->components_bags.size(); j++) {
        Table* final_table = calculate_tables(graph, decomp->components_bags[j], 
                                              postorder[j], optional_verts, variant);
        dom_set = get_solution(final_table); 
        
        printf("\n");
        for(auto it=dom_set->begin(); it!=dom_set->end(); it++) printf("soln set verts=%d\n", *it);
        
        component_domsets.push_back(dom_set); 
    }
    return component_domsets;
}


Table* calculate_tables(Graph* graph, std::vector<Set*>& bags, 
                        std::vector<po_bag>& postorder, 
                        Set* optional_verts, Variant variant) {
    /*     
     * Dynamic programming algorithm, for dominating set on graphs 
     * w/ bounded treewidth.
     * 
     * NOTE For NICE tree decompositions. For both annotated and regular dominating set. 
     */
    
    std::vector<Table*> tables;
    
    for(int i=0; i<postorder.size(); i++) {     // O(n -- # of bags)
        int bag_index = postorder[i].bag_index;
        int num_children = postorder[i].num_children;
        int parent_bag_index = postorder[i].parent_bag_index;
        
        printf("\n\n_____________________bag_index=%d, parent_bag_index=%d, num_children=%d \n", 
               bag_index, parent_bag_index, num_children);
        printf("Table size=%ld\n", tables.size());
        
        if(num_children==2) {          //-----------------------JOIN bag
            printf("JOIN BAG %d\n", bag_index);
            
            //update table here
            Table* left_child_table = tables[tables.size()-2];  
            Table* right_child_table = tables[tables.size()-1];
            
            update_join_table(right_child_table, left_child_table, 
                              optional_verts, variant);           // reuses the right child table
            
            tables.erase(tables.begin()+tables.size()-2);   // deletes the left child table. 

        } else if(num_children==1) {     //-------------------either INTRODUCE or FORGET bag
            int child_bag_index = postorder[i-1].bag_index;
            int child_bag_table_index = tables.size()-1;
            
            Set* parent_bag = bags[bag_index];
            Set* child_bag = bags[child_bag_index];
            
            string type;
            int v;
            
            if(parent_bag->size() > child_bag->size()) {        // introduce node
                printf("INTRO bag %d\n", bag_index);
                printf("Child bag index=%d\n", child_bag_index);
                for(auto it=parent_bag->begin(); it!=parent_bag->end(); it++) {
                    int u = *it;
                    if(!child_bag->contains(u)) v = u;
                }
                printf("intro vert=%d\n", v);
                
                Table* child_table = tables[child_bag_table_index];
                update_introduce_table(graph, child_table, child_bag, optional_verts, v, variant);
                
            } else if(parent_bag->size() < child_bag->size()) { // forget node
                printf("FORGET BAG %d\n", bag_index);
                printf("Child bag index=%d\n", child_bag_index);
                printf("Child bag table index=%d\n", child_bag_table_index);
            
                for(auto it=child_bag->begin(); it!=child_bag->end(); it++) {
                    int u = *it;
                    if(!parent_bag->contains(u)) v = u;
                }
                printf("forget vert=%d\n", v);
                
                //update table here
                Table* child_table = tables[child_bag_table_index];
                update_forget_table(child_table, optional_verts, v, variant);
                
            } else {
                printf("ERROR: in find_bagtype(), parent and child should not have same size.\n");
            }
            
        } else if(num_children==0){    //----------------------------LEAF bag
            //leaf bag, need to initialize table.
            printf("LEAF BAG %d\n", bag_index);
            
            Table* table = initialize_leaf_table(graph, bags[bag_index], optional_verts, variant);
            tables.push_back(table);
            
            std::string type = "Leaf";
            print_table(table, type);
            table->print_tablelookups();
            
        } else {
            printf("ERROR in number of children in nice decomp.\n");
        }
    }

    printf("\n\n\n__________________________TABLES_________________________\n");
    print_tables(tables);
    
    return tables[tables.size()-1]; //return the last table.
}


//---


Table* initialize_leaf_table(Graph* graph, Set* bag, Set* optional_verts, Variant variant) {
    /*
     * ni*3^ni
     * 
     * ni->(loop over all v in bag)*
     * 3^ni->(loop over at most 3^ni rows in the table)
     * 
     * Leaf tables can have more than 1 vertex.
     */
    Table* table = new Table();
    
    Row* r1;
    Row* r2;
    Row* r3;
        
    int i=0;
    for(auto it=bag->begin(); it!=bag->end(); it++) {  //3^ni time
        int v = *it;
        table->vertices.push_back(v);

        if(i==0) {
            //intitial rows for first vert.
            r1 = new Row();
            r1->append_coloring(IN_DOMSET);     //appends the coloring to the row's key.
            table->insert_row(r1);
            if(r1->coloring.size()==bag->size()) {  //if bag size==1
                r1->update_Ac(locally_valid_coloring(graph, optional_verts, 
                                                 r1, table->vertices, variant));                   //k^2
            }
                    
            r2 = new Row();
            r2->append_coloring(DOMINATED);
            table->insert_row(r2);
            if(r2->coloring.size()==bag->size()) {
                r2->update_Ac(locally_valid_coloring(graph, optional_verts, 
                                                 r2, table->vertices, variant));                   //k^2
            }
            
            r3 = new Row();
            r3->append_coloring(NOT_DOMINATED);
            table->insert_row(r3);
            if(r3->coloring.size()==bag->size()) {
                r3->update_Ac(locally_valid_coloring(graph, optional_verts, 
                                                 r3, table->vertices, variant));                   //k^2
            }
            
        } else {
            int tab_size = table->get_table_size();
            for(int j=0; j<tab_size; j++) {  //at most 3^ni
                //Row* r_update = table->table[j];    //no need to make actual copy
                Row* r_update = table->get_row(j);    //no need to make actual copy
                
                Row* r2 = new Row(r_update);
                r2->append_coloring(DOMINATED);
                table->insert_row(r2);
                if(r2->coloring.size()==bag->size()) {
                    r2->update_Ac(locally_valid_coloring(graph, optional_verts, 
                                                     r2, table->vertices, variant));               //k^2
                }
                
                Row* r3 = new Row(r_update);
                r3->append_coloring(NOT_DOMINATED);
                table->insert_row(r3);
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


void update_introduce_table(Graph* graph, Table* child_table, 
                            Set* child_bag, Set* optional_verts, 
                            int v, Variant variant) {
    /*
     * Updates the child table to be the parent. 
     */
    child_table->vertices.push_back(v);

    Set* neighbors_v = child_bag->set_intersection(graph->neighbors(v));  
    int table_size = child_table->get_table_size();
    
    for(int i=0; i<table_size; i++) {           //3^ni 
        Row* r_update = child_table->get_row(i);                  
                
        Row* r2 = new Row(r_update);                    //k, but also *n!! (when constructing soln currently)
        r2->append_coloring(DOMINATED);
        child_table->insert_row(r2);
        
        Row* r3 = new Row(r_update);                    //k
        r3->append_coloring(NOT_DOMINATED);
        child_table->insert_row(r3);
        
        //----r_update: x is IN_DOMSET=1
        int new_col_key = phi(r_update, neighbors_v, child_table->vertices, v);   //k
        printf("------------------new_col_key = %d\n", new_col_key);
        
        int A_phi;
        if(variant == Variant::Indep_Dom_Set) {             //Independent variant
            Set* independent_check = new Set();
            
            int incr = new_col_key;
            int index = 0;
            //TODO needs testing for longer new_col_keys 
            while (incr > 0) {                              //k
                //loops over single digits in int back to front.
                int c = incr%10;
                incr /= 10;
                printf("DIGIT=%d\n", c);
                
                if(c==IN_DOMSET) {
                    int size = child_table->vertices.size();
                    independent_check->insert(child_table->vertices[size-index]);
                }
                index++;
            }

            bool indep = check_independent(graph, independent_check); //k^2
            if(!indep) A_phi = INF;
            else A_phi = child_table->lookup_Ac(new_col_key) + 1;
        } else {
            A_phi = child_table->lookup_Ac(new_col_key) + 1;  //Dom Set variant
        }
        r_update->update_Ac(A_phi);
        child_table->update_row_add(r_update, IN_DOMSET);
        //-----
        
        //r2: x is DOMINATED=0
        bool is_xoptional = optional_verts->contains(v);  // in optional set or not?  
        
        //it's either not the annotated version or x is not optional
        if(!is_xoptional) {  
            //need to check that the introduced vertex w. the coloring of DOMINATED is justified.
            bool in_ds=false;
            for(auto it=neighbors_v->begin(); it!=neighbors_v->end(); it++) { //*k^2!!
                int index = child_table->get_vertex_col_index(*it);  
                
                //x has a neighbor IN_DOMSET
                if(r2->coloring[index]==IN_DOMSET) {
                    in_ds = true;
                }
            }
            if(!in_ds) r2->update_Ac(INF);
        } else { //annotated version and x IS an optional vertex, dont need to justify.
            //Ai(c x {DOMINATED}) <- Aj(c)  i.e. do nothing
        }
                
        //r3: x is NOT_DOMINATED=3  -- do nothing
    
    }
    delete neighbors_v;
    
    std::string type = "Introduce";
    print_table(child_table, type);
    child_table->print_tablelookups();
}


void update_forget_table(Table* child_table, Set* optional_verts, 
                         int v, Variant variant) {
    /*
     * Updates the child_table to be the parent. 
     * NOTE nothing (??) changes for the indpendent variant.
     */   
    int v_index = child_table->get_vertex_col_index(v);      //k
    int table_size = child_table->get_table_size();
    bool is_xoptional = optional_verts->contains(v);  // in optional set or not? 
    
    child_table->vertices.erase(child_table->vertices.begin()+v_index);   //delete v from verts vec. k
    
    int curr_index = 0;
    for(int j=0; j<table_size; j++) {                       //3^ni
        Row* row_child = child_table->get_row(curr_index);
        
        int color_v = row_child->coloring[v_index];
        int Aj = row_child->get_Ac();
        
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
                if(color_v==IN_DOMSET)  {
                    //Adds the forgotten vertex to the soln set if necessary
                    row_par->domset_verts->insert(v);
                }
                else if(color_v==DOMINATED) {
                    row_par->domset_verts->remove(v); 
                    row_par->domset_verts = row_par->domset_verts->set_union(row_child->domset_verts);
                }
                row_par->update_Ac(Aj);
            }    
        } else { //annotated version and x IS an optional vertex. 
            if(Aj <= Ai) {
                if(color_v==IN_DOMSET)  {
                    //Adds the forgotten vertex to the soln set if necessary
                    row_par->domset_verts->insert(v);
                } else {  // when color_v is DOMINATED or NOT_DOMINATED
                    row_par->domset_verts->remove(v); 
                    row_par->domset_verts = row_par->domset_verts->set_union(row_child->domset_verts);
                }
                row_par->update_Ac(Aj);
            }
        }
        
        curr_index++;
    }

    std::string type = "Forget";
    print_table(child_table, type);
    child_table->print_tablelookups();
}


void update_join_table(Table* rightchildk, Table* leftchildj, 
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
            
        minAi_c(rightchildk, leftchildj, optional_verts, row_k, row_j);        
    }
    
    std::string type = "Join";
    print_table(rightchildk, type);
    rightchildk->print_tablelookups();
}



//-- Dom set Helper functions


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
    Set* independent_check;
    
    int A_c = 0;
    bool valid = false;
    for(int i=0; i<vertices.size(); i++) {  //at most k
        int xt = vertices[i];
        int coloring_i = row->coloring[i];
        
        //If finding independent dominating set, create new set.
        if(variant == Variant::Indep_Dom_Set) independent_check = new Set();

        if(coloring_i == IN_DOMSET) {
            A_c++;
            
            if(variant == Variant::Indep_Dom_Set) independent_check->insert(xt);
        }
        
        if(coloring_i == DOMINATED) {  //if a vertex is set to DOMINATED
            bool actually_dominated = false;
            bool is_xtoptional = optional_verts->contains(xt);  // in optional set or not?
            
            //xt is not optional, or it's the annotated version and it is optional
            if(!is_xtoptional) {  
                for(int j=0; j<vertices.size(); j++) {  //at most k
                    int coloring_j = row->coloring[j];
                    
                    if(i!=j && coloring_j== IN_DOMSET) {
                        if(graph->adjacent(vertices[j], vertices[i])) actually_dominated = true;
                    }
                }
            } else actually_dominated = true;
            
            if(!actually_dominated) return INF;
        }
    }
    
    if(variant == Variant::Indep_Dom_Set) {
        bool indep = true;
        if(independent_check->size() >= 2) {
            indep = check_independent(graph, independent_check);
        }
        
        if(!indep) return INF;  // vertices are not independet
        delete independent_check;
    }
    
    
    return A_c;
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


void minAi_c(Table* childk, Table* childj, Set* optional_verts, Row* row_k, Row* row_j) {
    /*
     * NOTE: this is where storing the vertices vector is important. Look into not storing
     * it in table.
     *
     * NOTE: nothing (??) needs changing for independent variant. 
     */

    //First, find all possible c' and c'' which divide c.
    
    std::vector<int> c_prime;
    std::vector<int> c_primeprime;
    
    int num_ones = 0;   //number of 1's (IN_DOMSET's) in the coloring. 
    for(int i=0; i<row_k->coloring.size(); i++) { //k
        int c_t = row_k->coloring[i]; 
        //printf("c_t=%d\n", c_t);
        
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
                
        //printf("c'=%d, c''=%d\n", cprime_key, cprimeprime_key);
        int val;
        if(A_jcprime==INF || A_kprimeprime==INF) val = INF;
        else val = A_jcprime + A_kprimeprime - num_ones;
        
        if(val < min_Ai) {   
            min_Ai = val;
            A_jcprime_ind_final = A_jcprime_ind;
            A_kprimeprime_ind_final = A_kprimeprime_ind;
        }
    }
    
    //NOTE constructing solution here.
    Row* rj = childj->get_row(A_jcprime_ind_final);
    Row* rk = childk->get_row(A_kprimeprime_ind_final);
    
    row_k->domset_verts = rj->domset_verts->set_union(rk->domset_verts);
    row_k->update_Ac(min_Ai);
}


int phi(Row* row, Set* neighbors, std::vector<int> vertices, int introduced_v) {
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



