
#include "table.hpp"

#include "domset_exact.hpp"   //NOTE temporary for print


#include <iostream>
#include <fstream>
#include <limits>  //infinity


Table::Table(int lab) {
    /*
     * Table constructor. 
     * 
     * For NICE decomps
     */
    this->label = lab;
}


Table::~Table() {
    for(int i=0; i<table.size(); i++) delete table[i];
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


int Table::lookup(int key) {
    //Returns the INDEX of the row in the table.
    if(!table_lookups.contains(key)) return -1;
    return table_lookups[key];
}


void Table::insert_row(Row* r) {
    int index = table.size();
    table.push_back(r);  //adds row to the table vector.
        
    //key is the concatenated coloring, value is the index in table.
    table_lookups.insert(r->key, index);
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

