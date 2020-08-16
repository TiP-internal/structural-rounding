
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
            soln_row->domset_verts->insert(soln_row->coloring[i]);
        }
    }
    
    return soln_row->domset_verts;
}


std::vector<Set*> calc_domset(Graph* graph, TreeDecomp* decomp, bool annotated_version) {
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
        Table* final_table = calculate_tables(graph, decomp->components_bags[j], postorder[j]);
        dom_set = get_solution(final_table); 
        
        printf("\n");
        for(auto it=dom_set->begin(); it!=dom_set->end(); it++) printf("soln set verts=%d\n", *it);
        
        component_tables.push_back(dom_set); 
    }
    return component_tables;
}


Table* calculate_tables(Graph* graph, std::vector<Set*>& bags, std::vector<po_bag>& postorder) {
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
            
            //get child table indices. NOTE this could probably be improved 
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
            right_child_table->update_join_table(left_child_table, bag_index);
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
                child_table->update_introduce_table(graph, child_bag, v, bag_index);
            
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
                child_table->update_forget_table(v, bag_index);
                
            } else {
                printf("ERROR: in find_bagtype(), parent and child should not have same size.\n");
            }
            
        } else if(num_children==0){    //----------------------------LEAF bag
            //leaf bag, need to initialize table.
            printf("LEAF BAG %d\n", bag_index);
            
            Table* table = new Table(bag_index);
            table->initialize_leaf_table(graph, bags[bag_index]);
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
