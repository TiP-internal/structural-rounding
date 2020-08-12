
#include "table.hpp"

#include "domset_exact.hpp"   //NOTE temporary for print


#include <iostream>
#include <fstream>
#include <limits>  //infinity


Table::Table(Graph* graph, Set* bag, std::string type, int lab) {
    /*
     * Table constructor.
     */
    
    this->graph = graph;
    this->bag = bag;
    this->label = lab;
    this->table_type = type;
}


Table::~Table() {
    for(int i=0; i<table.size(); i++) delete table[i];
}


void Table::initialize_leaf_table() {
    /*
     * ni*3^ni
     * 
     * ni->(loop over all v in bag)*
     * 3^ni->(loop over at most 3^ni rows in the table)
     * 
     * Leaf tables can have more than 1 vertex.
     */
    
    Row* r1;
    Row* r2;
    Row* r3;
        
    int i=0;
    for(auto it=bag->begin(); it!=bag->end(); it++) {  //3^ni time
        int v = *it;
        vertices.push_back(v);
        printf("Vertex=%d\n", v);

        if(i==0) {
            //intitial rows for first vert.
            r1 = new Row();
            r1->append_coloring(IN_DOMSET);
            insert_row(r1);
            if(r1->coloring.size()==bag->size()) {  //if bag size==1
                r1->A_c = locally_valid_coloring(r1);                   //k^2
            }
                    
            r2 = new Row();
            r2->append_coloring(DOMINATED);
            insert_row(r2);
            if(r2->coloring.size()==bag->size()) {
                r2->A_c = locally_valid_coloring(r2);                   //k^2
            }
            
            r3 = new Row();
            r3->append_coloring(NOT_DOMINATED);
            insert_row(r3);
            if(r3->coloring.size()==bag->size()) {
                r3->A_c = locally_valid_coloring(r3);                   //k^2
            }
            
        } else {
            int tab_size = table.size();
            for(int j=0; j<tab_size; j++) {  //at most 3^ni
                Row* r_update = table[j];    //no need to make actual copy
                
                Row* r2 = new Row(r_update);
                r2->append_coloring(DOMINATED);
                insert_row(r2);
                if(r2->coloring.size()==bag->size()) {
                    r2->A_c = locally_valid_coloring(r2);               //k^2
                }
                
                Row* r3 = new Row(r_update);
                r3->append_coloring(NOT_DOMINATED);
                insert_row(r3);
                if(r3->coloring.size()==bag->size()) {
                    r3->A_c = locally_valid_coloring(r3);               //k^2
                }
                
                update_row_add(r_update, IN_DOMSET);
                if(r_update->coloring.size()==bag->size()) {
                    r_update->A_c = locally_valid_coloring(r_update);   //k^2
                }
            }
        }
        i++;
    }
    
    if(table.size()!=pow(3, vertices.size())) {
        printf("ERROR: in creating all colorings. table_size=%ld \
        3^vertices.size()=%d", table.size(), pow(3, vertices.size()));
    }
}


void Table::update_join_table(Table* leftchildj, Table* rightchildk) {
    /*
     * j is left
     * k is right
     * 
     * i is parent table
     */
    
    int tab_size = leftchildj->table.size();
    this->vertices = leftchildj->vertices;  //O(ni)
    
    this->left_child_table = leftchildj;
    this->right_child_table = rightchildk;
    
    for(int i=0; i<tab_size; i++) {  //4^ni
        Row* row_j = leftchildj->table[i];
        Row* row_k = rightchildk->table[i];
        
        Row* row_i = new Row();  //no copying here.
            
        int Aic = minAi_c(leftchildj, rightchildk, row_i, row_k, row_j);        
        row_i->A_c = Aic;
        insert_row(row_i);
    }
    
    print_table(this, 1);
    print_lookups(this);
}


void Table::update_forget_table(Table* child, int v) {
    /*
     * 
     */    
    this->vertices = child->vertices;           //O(ni)
    this->intro_or_forget_vertex=v;
    this->child_table = child;
    
    int v_index = get_vertex_col_index(v);  
    vertices.erase(vertices.begin()+v_index);  //delete v from verts vec. O(ni)
    
    for(int j=0; j<child->table.size(); j++) {  //3^ni
        //TODO only make copy if rowj/v not in table already.
        Row* rowj = new Row(child->table[j]);  //copy O(ni)
        
        int color_v = rowj->coloring[v_index];
        int Aj_prevcolor_ind = child->lookup(rowj->key);
        int Aj_prevcolor = child->table[Aj_prevcolor_ind]->A_c;
        
        update_row_delete(rowj, v_index);  //O(ni)
        
        //copying child table to parent table so must insert row. 
        if(!table_lookups.contains(rowj->key)) {
            insert_row(rowj);
        }
        
        //get the actual row in the parent table.
        int par_rowj_ind = lookup(rowj->key);
        Row* rowj_par = table[par_rowj_ind];  //not copying
        
        int Ai_currcolor_ind = lookup(rowj_par->key);
        int Ai = table[Ai_currcolor_ind]->A_c;
        
        if(color_v!=NOT_DOMINATED && Aj_prevcolor < Ai) {
            rowj_par->A_c = Aj_prevcolor;
        }
    }

    print_table(this, 2);
    print_lookups(this);
}

        
void Table::update_introduce_table(Table* child, int v) {
    /*
     * 
     */
    Set* neighbors_v = bag->set_intersection(graph->neighbors(v));
    int tab_size = child->table_lookups.size();
    
    this->vertices = child->vertices;               //O(ni)
    this->vertices.push_back(v);
    this->child_table = child;
    this->intro_or_forget_vertex = v;
    
    for(int i=0; i<tab_size; i++) {  //3^ni 
        Row* r_update = new Row(child->table[i]);   //O(ni)
        insert_row(r_update);
        
        Row* r2 = new Row(r_update);                //O(ni)
        r2->append_coloring(DOMINATED);
        insert_row(r2);
        
        Row* r3 = new Row(r_update);                //O(ni)
        r3->append_coloring(NOT_DOMINATED);
        insert_row(r3);
        
        //TODO Double check these.
        //x is IN_DOMSET
        int new_col_key = r_update->phi(neighbors_v->size()); //k
        int A_phi_ind = child->lookup(new_col_key);
        int A_phi = child->table[A_phi_ind]->A_c;
        
        r_update->A_c = A_phi+1;
        update_row_add(r_update, IN_DOMSET);

        //x is DOMINATED=0
        bool in_ds=false;
        for(auto it=neighbors_v->begin(); it!=neighbors_v->end(); it++) { //*k^2!!
            int index = get_vertex_col_index(*it); 
            
            //x has a neighbor IN_DOMSET
            if(r2->coloring[index]==IN_DOMSET) {
                in_ds = true;
            }
        }
        if(!in_ds) r2->A_c=INF;
                
        //r3: x is NOT_DOMINATED  -- do nothing
    
    }
    delete neighbors_v;
    print_table(this, 1);
    print_lookups(this);
}


//-----------------


int Table::minAi_c(Table* childj, Table* childk, Row* row_i, Row* row_k, Row* row_j) {
    /*
     * 
     */
    
    //First, find all possible c' and c'' which divide c.
    std::vector<std::string> c_prime;
    std::vector<std::string> c_primeprime;
    
    int num_ones = 0;   //number of 1's (IN_DOMSET's) in the coloring. 
    int k = 0; //number of 0's in the coloring.
    for(int i=0; i<row_j->coloring.size(); i++) { //ni
        
        int c_t = row_j->coloring[i]; 
        update_coloring_add(row_i, c_t);
        
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
            k++;
        
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
                    std::string cprimeprime_t_str = c_primeprime[t]; //O(ni)
                    
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
    for(int i=0; i<c_prime.size(); i++) {  
        int A_jcprime_ind = childj->lookup(stoi(c_prime[i]));
        int A_jcprime = childj->table[A_jcprime_ind]->A_c;
        
        int A_kprimeprime_ind = childk->lookup(stoi(c_primeprime[i]));
        int A_kprimeprime = childk->table[A_kprimeprime_ind]->A_c;
                
        int val;
        if(A_jcprime==INF || A_kprimeprime==INF) val = INF;
        else val = A_jcprime + A_kprimeprime - num_ones;
        
        if(val < min_Ai) min_Ai = val;
    }
    
    return min_Ai;
}


int Table::get_vertex_col_index(int v) {
    /*
     * TODO k time
     */
    int index = -1;
    for(int i=0; i<this->vertices.size(); i++) {
        if(v==vertices[i]) return i;
    }
    return index;
}


int Table::locally_valid_coloring(Row* row) {
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
    printf("locally valid coloring: ");
    int A_c = 0;
    bool valid = false;
    for(int i=0; i<vertices.size(); i++) {  //at most k
        int coloring_i = row->coloring[i];
        printf(" %d, ", coloring_i);

        if(coloring_i == IN_DOMSET) {
            A_c++;
        }
        
        if(coloring_i == DOMINATED) {  //if a vertex is set to DOMINATED
            bool actually_dominated = false;
            
            for(int j=0; j<vertices.size(); j++) {  //at most k
                int coloring_j = row->coloring[j];
                
                //another vert in bag is in domset and is neighbor.
                if(i != j && coloring_j== IN_DOMSET) {  
                    if(graph->adjacent(vertices[j], vertices[i])) actually_dominated = true;
                }
            }
            if(!actually_dominated) return INF;
        }
    }
    printf("\n\n");
    return A_c;
}


int Table::lookup(int key) {
    //Returns the INDEX of the row in the table.
    return table_lookups[key];
}


void Table::insert_row(Row* r) {
    int index = table.size();
    table.push_back(r);  //adds row to the table vector.
        
    //key is the concatenated coloring, value is the index in table.
    table_lookups.insert(r->key, index);
}


void Table::update_coloring_add(Row* r, int value) {
    //updates the key in the row, along with the coloring vect. 
    //does not add to the table/table_lookups 
    r->append_coloring(value);
}

void Table::update_row_add(Row* r, int value) {
    /* must update key in table_lookups whenever we add values
     * to the current coloring.
     * 
     * used for updating a coloring which is already in the table lookups. 
     */
    
    //1. erase
    int index = lookup(r->key);
    table_lookups.erase(r->key);
    
    //2. re-insert into table_lookups
    r->append_coloring(value);
    table_lookups.insert(r->key, index);  //same index, updated key
}


void Table::update_row_delete(Row* r, int v_index) {
    //1. erase
    int index = lookup(r->key);
    table_lookups.erase(r->key);
    r->remove_from_coloring(v_index);  //O(ni)
}

