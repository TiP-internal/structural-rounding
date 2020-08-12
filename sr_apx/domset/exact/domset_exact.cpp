
#include "domset_exact.hpp"

#include <cstdio>

//---for testing
bool is_domset(Graph* graph, std::vector<int> domset) {
    for(auto it=graph->begin(); it!=graph->end(); it++) {
        int v=*it;
        bool adjacent = false;
        for(int i=0; i<domset.size(); i++) {
            int u=domset[i];
            if(graph->adjacent(v, u)) {
                adjacent=true;
            }
        }
        if(!adjacent) return false;
    }
    return true;
}

void print_table(Table* tab, int i) {
    printf("\n===========================================\n");
    printf("IN_DOMSET=1=1,   DOMINATED=2=0,   NOT_DOMINATED=3=0hat\n");
    printf("          Table: %d, Label: %d, Type: %s  \n", i, tab->label, tab->table_type.c_str());
    printf("  vertices \t\t |  A_ci \t | Soln. Set \n");
    printf("|");
    for(int j=0; j<tab->vertices.size(); j++) {
        printf("  %d  |", tab->vertices[j]);
    }
    printf("\n___________________________________________\n");
    printf("\n");
    
    for(int j=0; j<tab->table.size(); j++) {
        Row* r = tab->table[j];
        printf("|");
        for(int k=0; k<r->coloring.size(); k++) {
            printf("  %d  |", r->coloring[k]);
        }
        printf("%15d \t | ", r->A_c);
        for(auto it=r->domset_verts->begin(); it!=r->domset_verts->end(); it++) {
            printf(" %d,", *it);
        }
        printf("\n");
    }
    
    printf("\n===========================================\n\n\n");
}

void print_tables(std::vector<Table*> tables) {
    for(int i=0; i<tables.size(); i++) {
        Table* tab = tables[i];
        print_table(tab, i);
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
     * Traverses through tables from root table to leafs to get solution.
     * TODO
     */
    Set* dom_set = new Set();
    //calc answer

    return dom_set;
}


Set* calc_domset(Graph* graph, TreeDecomp* decomp) {
    /*
     * Calculates the minimum dominating set given a tree decomp.
     * Must calculate for each component in the decomp. 
     * 
     * TODO must update to handle each component.
     */
    std::vector<std::vector<po_bag>> postorder = decomp->get_post_order();

    Set* dom_set;
    // std::vector<Table*> component_tables;
    for(int j=0; j<decomp->components_bags.size(); j++) {
        Table* final_table = calculate_tables(graph, decomp->components_bags[j], postorder[j]);
        dom_set = get_solution(final_table);    
        
        // component_tables.push_back(comp); 
    }
    return dom_set;
}


Table* calculate_tables(Graph* graph, std::vector<Set*>& bags, std::vector<po_bag>& postorder) {
    /*     
     * Dynamic programming algorithm, for dominating set on graphs 
     * w/ bounded treewidth.
     * 
     * NOTE For NICE tree decompositions. 
     */
    std::vector<Table*> tables;
    
    for(int i=0; i<postorder.size(); i++) { // O(n)
        int bag_index = postorder[i].bag_index;
        int num_children = postorder[i].num_children;
        int parent_bag_index = postorder[i].parent_bag_index;
        
        printf("\n\n_____________________bag_index=%d, parent_bag_index=%d, num_children=%d \n", 
               bag_index, parent_bag_index, num_children);
        
        if(num_children==2) {          //-----------------------JOIN bag
            printf("JOIN BAG %d\n", bag_index);
            
            Table* table = new Table(graph, bags[bag_index], "join_bag", bag_index);
            tables.push_back(table);
            
            //get child table indices. TODO this could probably be improved 
            int table_index_child_left=-99;
            int table_index_child_right=-99;
            for(int j=0; j<postorder.size(); j++) {             //O(n)
                if(postorder[j].parent_bag_index==bag_index) {
                    if(table_index_child_left==-99) {
                        table_index_child_left=j;
                    } else if(table_index_child_right==-99) {
                        table_index_child_right=j;
                    }
                }
            }
            
            printf("left child ind=%d, right child ind=%d\n", 
                   table_index_child_left, table_index_child_right);
            
            //update table here
            table->update_join_table(tables[table_index_child_left],
                                     tables[table_index_child_right]);
            
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
                
                Table* table = new Table(graph, parent_bag, "introduce_bag", bag_index);
                tables.push_back(table);
                
                table->update_introduce_table(tables[child_bag_table_index], v);
            
            } else if(parent_bag->size() < child_bag->size()) { // forget node
                printf("FORGET BAG %d\n", bag_index);
                printf("Child bag index=%d\n", child_bag_index);
            
                for(auto it=child_bag->begin(); it!=child_bag->end(); it++) {
                    int u = *it;
                    if(!parent_bag->contains(u)) v = u;
                }
                
                Table* table = new Table(graph, parent_bag, "forget_bag", bag_index);
                tables.push_back(table);
                
                //update table here
                table->update_forget_table(tables[child_bag_table_index], v);
                
            } else {
                printf("ERROR: in find_bagtype(), parent and child should not have same size.\n");
            }
            
            
        } else if(num_children==0){    //----------------------------LEAF bag
            //leaf bag, need to initialize table.
            printf("LEAF BAG %d\n", bag_index);
            
            Table* table = new Table(graph, bags[bag_index], "leaf_bag", bag_index);
            table->initialize_leaf_table();
            tables.push_back(table);
            
            print_table(table, 0);
            print_lookups(table);
            
        } else {
            printf("ERROR in number of children in nice decomp.\n");
        }
    }

    //printf("\n\n\n__________________________TABLES_________________________\n");
    //print_tables(tables);
    
    return tables[tables.size()-1]; //return the last table.
}
