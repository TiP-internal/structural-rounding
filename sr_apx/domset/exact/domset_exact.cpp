
#include "domset_exact.hpp"

#include <cstdio>

//---for testing
bool is_domset(Graph* graph, Set* domset) {
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int v=*it;
        bool adjacent = false;
        //for(int i=0; i<domset->size(); i++) {
        for(auto ids=domset->begin(); ids!=domset->end(); ids++) {
            int u=*ids;
            if(graph->adjacent(v, u)) {
                adjacent=true;
            }
        }
        if(!adjacent) return false;
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
        if(!adjacent) return false;
    }
    return true;
}

void print_row(Row* row) {
    printf("|");
    for(int k=0; k<row->coloring.size(); k++) {
        printf("  %d  |", row->coloring[k]);
    }
    printf("%15d \t | ", row->A_c);
    for(auto it=row->domset_verts->begin(); it!=row->domset_verts->end(); it++) {
        printf(" %d,", *it);
    }
    printf("\n");
}

void print_table(Table* tab, std::string table_type) {
    printf("\n===========================================\n");
    printf("IN_DOMSET=1=1,   DOMINATED=2=0,   NOT_DOMINATED=3=0hat\n");
    printf("            Label: %d, Type: %s  \n", tab->label, table_type.c_str());
    printf("  vertices \t\t |  A_ci \t | Soln. Set \n");
    printf("|");
    for(int j=0; j<tab->vertices.size(); j++) {
        printf("  %d  |", tab->vertices[j]);
    }
    printf("\n___________________________________________\n");
    printf("\n");
    
    for(int j=0; j<tab->table.size(); j++) {
        Row* r = tab->table[j];
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

void print_lookups(Table* table) {
    printf("\n===============Table Lookups==================\n");
    for(auto it=table->table_lookups.begin(); it!=table->table_lookups.end(); it++) {
        int key = *it;
        printf("coloring=%d, \t table_index=%d\n", key, table->table_lookups[key]);
    }
    
    printf("\n===========================================\n\n\n");
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
    for(int i=0; i<table->table.size(); i++) {
        Row* row = table->table[i];
        printf("row->Ai=%d\n", row->A_c);
        
        if(row->A_c < min_Ai) {
            min_Ai=row->A_c;
            soln_index=i;
        }
    }
    
    if(soln_index==-1) printf("ERROR: soln row not found.\n");
        
    //Add verts (set as IN_DOMSET) in the coloring to the soln set. 
    Row* soln_row = table->table[soln_index];
    for(int i=0; i<soln_row->coloring.size(); i++) {
        if(soln_row->coloring[i]==IN_DOMSET) {
            soln_row->domset_verts->insert(table->vertices[i]);
        }
    }
    
    return soln_row->domset_verts;
}


std::vector<Set*> calc_domset(Graph* graph, TreeDecomp* decomp, 
                              Set* optional_verts, bool annotated_version) {
    /*
     * Calculates the minimum dominating set given a tree decomp.
     * Must calculate for each component in the decomp. 
     * 
     * Calculates annotated and regular version (separately) of the Dominating
     * Set problem. 
    */
    std::vector<std::vector<po_bag>> postorder = decomp->get_post_order();

    Set* dom_set;
    std::vector<Set*> component_tables;
    for(int j=0; j<decomp->components_bags.size(); j++) {
        Table* final_table = calculate_tables(graph, decomp->components_bags[j], 
                                              postorder[j], optional_verts,
                                              annotated_version);
        dom_set = get_solution(final_table); 
        
        printf("\n");
        for(auto it=dom_set->begin(); it!=dom_set->end(); it++) printf("soln set verts=%d\n", *it);
        
        component_tables.push_back(dom_set); 
    }
    return component_tables;
}


Table* calculate_tables(Graph* graph, std::vector<Set*>& bags, 
                        std::vector<po_bag>& postorder, 
                        Set* optional_verts, bool annotated_version) {
    /*     
     * Dynamic programming algorithm, for dominating set on graphs 
     * w/ bounded treewidth.
     * 
     * NOTE For NICE tree decompositions. For both annotated and regular dominating set. 
     */
    std::vector<Table*> tables;
    
    for(int i=0; i<postorder.size(); i++) { // O(n)
        int bag_index = postorder[i].bag_index;
        int num_children = postorder[i].num_children;
        int parent_bag_index = postorder[i].parent_bag_index;
        
        printf("\n\n_____________________bag_index=%d, parent_bag_index=%d, num_children=%d \n", 
               bag_index, parent_bag_index, num_children);
        printf("Table size=%ld\n", tables.size());
        
        if(num_children==2) {          //-----------------------JOIN bag
            printf("JOIN BAG %d\n", bag_index);
            
            //get child table indices, this could probably be improved 
            int table_index_child_left=-99;
            int table_index_child_right=-99;
//             for(int j=0; j<postorder.size(); j++) {             //O(n)
//                 if(postorder[j].parent_bag_index==bag_index) {
//                     if(table_index_child_left==-99) {
//                         table_index_child_left=j;
//                     } else if(table_index_child_right==-99) {
//                         table_index_child_right=j;
//                     }
//                 }
//             }
            
//             printf("left child ind=%d, right child ind=%d\n", 
//                    table_index_child_left, table_index_child_right);
            
            //update table here
            Table* left_child_table = tables[tables.size()-2];  //TODO find correct left table
            Table* right_child_table = tables[tables.size()-1];
            
            update_join_table(right_child_table, left_child_table, 
                              optional_verts, bag_index, 
                              annotated_version);
            
            tables.erase(tables.begin()+tables.size()-2);  // delete the left child table. 

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
                
                Table* child_table = tables[child_bag_table_index];
                update_introduce_table(graph, child_table, child_bag, optional_verts,
                                       v, bag_index, annotated_version);
                
            } else if(parent_bag->size() < child_bag->size()) { // forget node
                printf("FORGET BAG %d\n", bag_index);
                printf("Child bag index=%d\n", child_bag_index);
                printf("Child bag table index=%d\n", child_bag_table_index);
            
                for(auto it=child_bag->begin(); it!=child_bag->end(); it++) {
                    int u = *it;
                    if(!parent_bag->contains(u)) v = u;
                }
                
                //update table here
                Table* child_table = tables[child_bag_table_index];
                update_forget_table(child_table, optional_verts, 
                                    v, bag_index, annotated_version);
                
            } else {
                printf("ERROR: in find_bagtype(), parent and child should not have same size.\n");
            }
            
        } else if(num_children==0){    //----------------------------LEAF bag
            //leaf bag, need to initialize table.
            printf("LEAF BAG %d\n", bag_index);
            
            Table* table = initialize_leaf_table(graph, bags[bag_index], optional_verts,
                                                 bag_index, annotated_version);
            tables.push_back(table);
            
            std::string type = "Leaf";
            print_table(table, type);
            print_lookups(table);
            
        } else {
            printf("ERROR in number of children in nice decomp.\n");
        }
    }

    printf("\n\n\n__________________________TABLES_________________________\n");
    print_tables(tables);
    
    return tables[tables.size()-1]; //return the last table.
}


//---


Table* initialize_leaf_table(Graph* graph, Set* bag, Set* optional_verts,
                             int bag_index, bool annotated_version) {
    /*
     * ni*3^ni
     * 
     * ni->(loop over all v in bag)*
     * 3^ni->(loop over at most 3^ni rows in the table)
     * 
     * Leaf tables can have more than 1 vertex.
     */
    Table* table = new Table(bag_index);
    
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
            r1->append_coloring(IN_DOMSET);
            table->insert_row(r1);
            if(r1->coloring.size()==bag->size()) {  //if bag size==1
                r1->A_c = locally_valid_coloring(graph, optional_verts, 
                                                 r1, table->vertices, 
                                                 annotated_version);                   //k^2
            }
                    
            r2 = new Row();
            r2->append_coloring(DOMINATED);
            table->insert_row(r2);
            if(r2->coloring.size()==bag->size()) {
                r2->A_c = locally_valid_coloring(graph, optional_verts, 
                                                 r2, table->vertices, 
                                                 annotated_version);                   //k^2
            }
            
            r3 = new Row();
            r3->append_coloring(NOT_DOMINATED);
            table->insert_row(r3);
            if(r3->coloring.size()==bag->size()) {
                r3->A_c = locally_valid_coloring(graph, optional_verts, 
                                                 r3, table->vertices, 
                                                 annotated_version);                   //k^2
            }
            
        } else {
            int tab_size = table->table.size();
            for(int j=0; j<tab_size; j++) {  //at most 3^ni
                Row* r_update = table->table[j];    //no need to make actual copy
                
                Row* r2 = new Row(r_update);
                r2->append_coloring(DOMINATED);
                table->insert_row(r2);
                if(r2->coloring.size()==bag->size()) {
                    r2->A_c = locally_valid_coloring(graph, optional_verts, 
                                                     r2, table->vertices, 
                                                     annotated_version);               //k^2
                }
                
                Row* r3 = new Row(r_update);
                r3->append_coloring(NOT_DOMINATED);
                table->insert_row(r3);
                if(r3->coloring.size()==bag->size()) {
                    r3->A_c = locally_valid_coloring(graph, optional_verts, 
                                                     r3, table->vertices, 
                                                     annotated_version);               //k^2
                }
                
                table->update_row_add(r_update, IN_DOMSET);
                if(r_update->coloring.size()==bag->size()) {
                    r_update->A_c = locally_valid_coloring(graph, optional_verts, 
                                                           r_update, table->vertices, 
                                                           annotated_version);   //k^2
                }
            }
        }
        i++;
    }
    
    if(table->table.size()!=pow(3, table->vertices.size())) {
        printf("ERROR: in creating all colorings. table_size=%ld \
        3^vertices.size()=%d", table->table.size(), pow(3, table->vertices.size()));
    }
    
    return table;
}


void update_introduce_table(Graph* graph, Table* child_table, Set* child_bag, Set* optional_verts,
                            int v, int new_label, bool annotated_version) {
    /*
     * Updates the child table to be the parent. 
     */
    child_table->label = new_label;
    child_table->vertices.push_back(v);

    //NOTE could just add bag as param and make sure to insert the child's bag into fun call.    
    Set* neighbors_v = child_bag->set_intersection(graph->neighbors(v));  
    int tab_size = child_table->table_lookups.size();
    
    for(int i=0; i<tab_size; i++) {  //3^ni 
        Row* r_update = child_table->table[i];                   //k
                
        Row* r2 = new Row(r_update);                //k
        r2->append_coloring(DOMINATED);
        child_table->insert_row(r2);
        
        Row* r3 = new Row(r_update);                //k
        r3->append_coloring(NOT_DOMINATED);
        child_table->insert_row(r3);
        
        //x is IN_DOMSET
        //int new_col_key = r_update->phi(neighbors_v->size()); //k
        int new_col_key = phi(r_update, neighbors_v->size()); //k
        int A_phi_ind = child_table->lookup(new_col_key);
        int A_phi = child_table->table[A_phi_ind]->A_c;
        
        r_update->A_c = A_phi+1;
        child_table->update_row_add(r_update, IN_DOMSET);

        //x is DOMINATED=0
        bool in_ds=false;
        for(auto it=neighbors_v->begin(); it!=neighbors_v->end(); it++) { //*k^2!!
            int index = child_table->get_vertex_col_index(*it); 
            
            //x has a neighbor IN_DOMSET
            if(r2->coloring[index]==IN_DOMSET) {
                in_ds = true;
            }
        }
        if(!in_ds) r2->A_c=INF;
                
        //r3: x is NOT_DOMINATED  -- do nothing
    
    }
    delete neighbors_v;
    
    std::string type = "Introduce";
    print_table(child_table, type);
    print_lookups(child_table);
}


void update_forget_table(Table* child_table, Set* optional_verts, 
                         int v, int new_label, bool annotated_version) {
    /*
     * Updates the child_table to be the parent. 
     */   
    child_table->label = new_label;

    //NOTE this extra loop could potentially be removed if we could
    // guarantee that the "forgotten" vertex is at the first index.
    int v_index = child_table->get_vertex_col_index(v);      //k
    child_table->vertices.erase(child_table->vertices.begin()+v_index);   //delete v from verts vec. k
    int table_size = child_table->table.size();
    
    int curr_index = 0;
    for(int j=0; j<table_size; j++) {                   //3^ni
        Row* row_child = child_table->table[curr_index];
        
        int color_v = row_child->coloring[v_index];
        int Aj = row_child->A_c;
        
        child_table->table_lookups.erase(row_child->key);
        row_child->remove_from_coloring(v_index);       //k

        int par_row_ind;
        Row* row_par;
        if(child_table->table_lookups.contains(row_child->key)) {  // row has already been added
            //delete_row(curr_index);
            child_table->table.erase(child_table->table.begin()+curr_index); 
            curr_index--;
            
            par_row_ind = child_table->lookup(row_child->key);   // get the parent row which has been added.
            row_par = child_table->table[par_row_ind];
        } else {  //row needs readding to table lookups
            par_row_ind = curr_index;
            row_par = child_table->table[par_row_ind];  
            
            //add to table lookups.
            child_table->table_lookups.insert(row_par->key, par_row_ind);  //same index, updated key
        }
        
        int Ai_ind = child_table->lookup(row_par->key);
        int Ai = child_table->table[Ai_ind]->A_c;
        
        //NOTE may need to add extra constraint here
        if(color_v!=NOT_DOMINATED && Aj <= Ai) {
            if(color_v==IN_DOMSET)  {
                //Adds the forgotten vertex to the soln set if necessary
                row_par->domset_verts->insert(v);
            }
            else if(color_v==DOMINATED) {
                row_par->domset_verts->remove(v); 
                row_par->domset_verts = row_par->domset_verts->set_union(row_child->domset_verts);
            }
            
            row_par->A_c = Aj;
        }
        curr_index++;
    }

    std::string type = "Forget";
    print_table(child_table, type);
    print_lookups(child_table);
}


void update_join_table(Table* rightchildk, Table* leftchildj, Set* optional_verts,
                       int new_label, bool annotated_version) {
    /*
     * j is left  
     * k is right  -- calling table
     * 
     * i is parent table
     * 
     * Updates the rightchild table to be the parent. 
     */
    rightchildk->label = new_label;
    int tab_size = rightchildk->table.size();
    
    for(int i=0; i<tab_size; i++) {  //4^ni
        Row* row_j = leftchildj->table[i];
        Row* row_k = rightchildk->table[i];
            
        minAi_c(rightchildk, leftchildj, row_k, row_j);        
    }
    
    std::string type = "Join";
    print_table(rightchildk, type);
    print_lookups(rightchildk);
}



//-- Dom set Helper functions


int locally_valid_coloring(Graph* graph, Set* optional_verts, Row* row, 
                           std::vector<int> &vertices, bool annotated_version) {
    /*
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
    int A_c = 0;
    bool valid = false;
    for(int i=0; i<vertices.size(); i++) {  //at most k
        int coloring_i = row->coloring[i];

        if(coloring_i == IN_DOMSET) {
            A_c++;
        }
        
        if(coloring_i == DOMINATED) {  //if a vertex is set to DOMINATED
            bool actually_dominated = false;
                
            for(int j=0; j<vertices.size(); j++) {  //at most k
                int coloring_j = row->coloring[j];
                
                if(i!=j) {
                    //another vert in bag is in domset and is neighbor.
                    if(!annotated_version) {
                        if(coloring_j== IN_DOMSET) {  
                            if(graph->adjacent(vertices[j], vertices[i])) actually_dominated = true;
                        }
                    } else {  //NOTE Changed for annotated version  TODO double check this
                        if(optional_verts->contains(vertices[i])) { //CAN be dominated
                            if(graph->adjacent(vertices[j], vertices[i])) actually_dominated = true;
                        } else {  //not in opt verts set. MUST be dominated
                            if(coloring_j== IN_DOMSET) {
                                if(graph->adjacent(vertices[j], vertices[i])) actually_dominated = true;
                            }
                        }
                    }
                }
            }
            
            if(!actually_dominated) return INF;
        }
    }
    return A_c;
}


void minAi_c(Table* childk, Table* childj, Row* row_k, Row* row_j) {
    /*
     * j is this child
     */

    //First, find all possible c' and c'' which divide c.
    std::vector<std::string> c_prime;
    std::vector<std::string> c_primeprime;
    
    int num_ones = 0;   //number of 1's (IN_DOMSET's) in the coloring. 
    for(int i=0; i<row_k->coloring.size(); i++) { //k
        int c_t = row_k->coloring[i]; 
        //printf("c_t=%d\n", c_t);
        
        if(c_t == IN_DOMSET) num_ones++;
        
        int cprime_size = c_prime.size();
        int cprimeprime_size = c_primeprime.size();
        
        if(cprime_size != cprimeprime_size) printf("ERROR: in divide().\n");
        
        if(c_t == IN_DOMSET || c_t==NOT_DOMINATED) {
            if(cprime_size==0 && cprimeprime_size==0) {
                std::string cprime0_str = "";
                cprime0_str = cprime0_str+std::to_string(c_t);
                c_prime.push_back(cprime0_str);    //c_t' == c_t
                
                std::string cprimeprime0_str = "";
                cprimeprime0_str = cprimeprime0_str+std::to_string(c_t);
                c_primeprime.push_back(cprimeprime0_str);   //c_t'' == c_t
            } else {
                for(int t=0; t<cprime_size; t++) {
                    //c_t' == c_t
                    c_prime[t] = c_prime[t] + std::to_string(c_t);   
                    
                    //c_t'' == c_t
                    c_primeprime[t] = c_primeprime[t] + std::to_string(c_t); 
                }
            }
        } else if (c_t == DOMINATED ) {
            if(cprime_size==0 && cprimeprime_size==0) {
                //---- c'_t = 0
                std::string cprime0_str = "";
                cprime0_str = cprime0_str+std::to_string(DOMINATED);
                c_prime.push_back(cprime0_str);
                
                //c''_t = 0hat
                std::string cprimeprime0_str = "";
                cprimeprime0_str = cprimeprime0_str+std::to_string(NOT_DOMINATED);
                c_primeprime.push_back(cprimeprime0_str);
                
                //---- c'_t = 0hat
                std::string cprime1_str = "";
                cprime1_str = cprime1_str+std::to_string(NOT_DOMINATED);
                c_prime.push_back(cprime1_str);
                
                //c''_t = 0
                std::string cprimeprime1_str = "";
                cprimeprime1_str = cprimeprime1_str+std::to_string(DOMINATED);
                c_primeprime.push_back(cprimeprime1_str);
            } else {
                //NOTE combine these for loops
                for(int t=0; t<cprime_size; t++) {
                    //1. make copy
                    std::string cprime_t_str = c_prime[t];
                    
                    //2. update current
                    c_prime[t] = c_prime[t]+std::to_string(DOMINATED);
                    
                    //3. add to second pair
                    cprime_t_str = cprime_t_str+std::to_string(NOT_DOMINATED);
                    
                    //4. add second pair to end of vector
                    c_prime.push_back(cprime_t_str);
                    
                    //---------
                    //1. make copy
                    std::string cprimeprime_t_str = c_primeprime[t]; //k
                    
                    //2. update current
                    c_primeprime[t]=c_primeprime[t]+std::to_string(NOT_DOMINATED);
                    
                    //3. add to second pair
                    cprimeprime_t_str = cprimeprime_t_str+std::to_string(DOMINATED);
                    
                    //4. add second pair to end of vector
                    c_primeprime.push_back(cprimeprime_t_str);
                }
            }
        }
    }
    
    //Now, find Ai(c) ← min{Aj(c')+ Ak(c'') − #1(c) | c' and c'' divide c }
    int min_Ai = INF;
    int A_jcprime_ind_final = -1;
    int A_kprimeprime_ind_final = -1;

    for(int i=0; i<c_prime.size(); i++) {  
        int A_jcprime_ind = childj->lookup(stoi(c_prime[i]));
        int A_jcprime = childj->table[A_jcprime_ind]->A_c;
        
        int A_kprimeprime_ind = childk->lookup(stoi(c_primeprime[i]));
        int A_kprimeprime = childk->table[A_kprimeprime_ind]->A_c;
                
        //printf("c'=%d, c''=%d\n", stoi(c_prime[i]), stoi(c_primeprime[i]));
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
    Row* rj = childj->table[A_jcprime_ind_final];
    Row* rk = childk->table[A_kprimeprime_ind_final];
    
    row_k->domset_verts = rj->domset_verts->set_union(rk->domset_verts);
    row_k->A_c = min_Ai;
}


int phi(Row* row, int num_neigbs) {
   /*
    *  φ : {0, ˆ0, 1}^nj → {0, ˆ0, 1}^nj 
    * on the set of colorings of Xj. 
    * For c =(c1,... ,c_nj ) ∈ {0, ˆ0, 1}^nj , let φ(c):= (c'_1,... ,c'_nj ) 
    * such that
    * 
    * c'_t = 0ˆ if t ∈ {p1,... ,ps} and ct =0  OR
    * c'_t = c_t otherwise.
    * 
    * TODO use actual indexes of neighbors in the coloring vertex?
    * 
    * returns key of the new coloring
    */
    
    std::string k_str = "";
    for(int i=0; i<row->coloring.size(); i++) {
        if( row->coloring[i] == DOMINATED  && i< num_neigbs) {
            std::string k = std::to_string(NOT_DOMINATED); 
            k_str = k_str+k;
        } else {
            std::string k = std::to_string(row->coloring[i]); 
            k_str = k_str+k;
        }
    }
    return stoi(k_str);
}



