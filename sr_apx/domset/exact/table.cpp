
#include "table.hpp"

#include "domset_exact.hpp"   //NOTE temporary for print


#include <iostream>
#include <fstream>
#include <limits>  //infinity


Row::Row() {
    this->A_c = INF;
    this->key = -1;
    
    this->childl_ind = -1;
    this->childr_ind = -1;
}

Row::Row(const Row* r) {  //copy constructor
    for(int i=0; i<r->coloring.size(); i++) {
        this->coloring.push_back(r->coloring[i]);
    }
    this->A_c = r->A_c;
    this->key = r->key;
    
    this->childl_ind=r->childl_ind;
    this->childr_ind=r->childr_ind;
}

Row::~Row() {}

int Row::get_Ac() {
    return A_c;
}

void Row::update_Ac(int newAc) {
    A_c = newAc;
}

int Row::get_key() {
    return key;
}

void Row::update_key(int newKey) {
    key = newKey;
}

void Row::append_coloring(int val) {
    //val should only ever be a single int (1,2, or 3)
    coloring.push_back(val);
    
    if (key==-1) key = val;
    else {
        key = key*10;  
        key = key+val;
    }
}

void Row::remove_from_coloring(int v_index) {
    //4*k time
    //This is only called once in the forget table.
    coloring.erase(coloring.begin()+v_index);  //erase the coloring for v
    
    std::string key_str = std::to_string(key);  //k
    key_str.erase(v_index, 1); //deletes one characters at index?
    
    //Complexity: "Unspecified, but generally linear in the number 
    //of characters interpreted."
    key = stoi(key_str);
}


/* The below functions are what track how the solutions are created.
 * ie are used during the second pass of the algorithm to construct the
 * solution set.
 */
int Row::get_childl_table_ind() {
    return childl_ind;
}

void Row::set_childl_table_ind(int childl_ind) {
    this->childl_ind=childl_ind;
}

int Row::get_childr_table_ind() {
    return childr_ind;
}

void Row::set_childr_table_ind(int childr_ind) {
    this->childr_ind=childr_ind;
}



Table::Table() {}

Table::~Table() {
    for(int i=0; i<table.size(); i++) {
        delete table[i];
    }
}


void Table::set_pobag(po_bag po) {
    tables_pobag=po;
}

po_bag Table::get_pobag() {
    return tables_pobag;
}

int Table::get_vertex_col_index(int v) {
    /*
     * k time
     */
    int index = -1;
    for(int i=0; i<this->vertices.size(); i++) {
        if(v==vertices[i]) return i;
    }
    return index;
}


int Table::lookup_table_index(int key) {
    //Returns the INDEX of the row in the table.
    if(!table_lookups.contains(key)) return -1;
    return table_lookups[key];
}


Row* Table::lookup_row(int key) {
    //Returns the row from the table w. the given key
    int index = lookup_table_index(key);
    return get_row(index);
}

Row* Table::get_row(int index) {
    //Returns the row from the table at the given index.
    if(index >= table.size() || index < 0) {
        throw("ERROR: out of bounds table index=%d, table size=%ld\n", index, table.size());
    }
    return table[index];
}


Row* Table::create_row(Row* r_update, int coloring) {
    Row* r = new Row(r_update);        //k
    if(coloring!=-1) r->append_coloring(coloring);
    insert_row(r);
    
    return r;
}

Row* Table::create_row(int coloring) {
    Row* r = new Row();        
    if(coloring!=-1) r->append_coloring(coloring);
    insert_row(r);
    
    return r;
}

int Table::lookup_Ac(int key) {
    /*
     * Given a key, returns the Ac value of the row corresponding
     * to that key.
     */
    int index = lookup_table_index(key);
    return table[index]->get_Ac();
}

void Table::insert_row(Row* r) {
    //only inserts if r not already inserted 
    if(lookup_table_index(r->get_key()) == -1) {
        int index = table.size();
        table.push_back(r);         //adds row to the table vector.
            
        //key is the concatenated coloring, value is the index in table.
        //table_lookups.insert(r->get_key(), index);
        table_lookups_insert(r->get_key(), index);
    }
}

void Table::delete_row(int index) {
    table.erase(table.begin()+index);  //*3^ni
}

void Table::table_lookups_insert(int key, int table_index) {
    if(key <0) throw("Incorrect key: table_lookups_insert\n");
    
    //Adds an entry to table_lookups_insert Map
    table_lookups.insert(key, table_index);
}


void Table::table_lookups_remove(int key) {
    //Removes an entry from the table_lookups Map
    table_lookups.erase(key);
}

bool Table::table_lookups_contains(int key) {
    return table_lookups.contains(key);
}


void Table::update_row_add(Row* r, int value) {
    /* must update key in table_lookups whenever we add values
     * to the current coloring.
     * 
     * used for updating a coloring which is already in the table lookups. 
     */
    //1. erase
    int index = lookup_table_index(r->get_key());
    table_lookups_remove(r->get_key());
    
    //2. re-insert into table_lookups
    r->append_coloring(value);
    //table_lookups.insert(r->get_key(), index);  //same index, updated key
    table_lookups_insert(r->get_key(), index);
}

int Table::get_table_size() {
    return table_lookups.size();
}


