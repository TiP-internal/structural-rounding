
#include "table.hpp"

#include "domset_exact.hpp"   //NOTE temporary for print


#include <iostream>
#include <fstream>
#include <limits>  //infinity


Row::Row() {
    A_c = INF;
    key = -1;
    
    childl_ind = -1;
    childr_ind = -1;
}

Row::Row(const Row &r) {  //copy constructor
    coloring = r.coloring;
    A_c = r.A_c;
    key = r.key;
    
    childl_ind=r.childl_ind;
    childr_ind=r.childr_ind;
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

void Row::set_childl_table_ind(int childl_ind_new) {
    childl_ind=childl_ind_new;
}

int Row::get_childr_table_ind() {
    return childr_ind;
}

void Row::set_childr_table_ind(int childr_indnew) {
    childr_ind=childr_indnew;
}




Table::Table() {}

Table::~Table() {}

void Table::set_pobag(po_bag po) {
    tables_pobag=po;
}

po_bag Table::get_pobag() {
    return tables_pobag;
}

void Table::update_Ac(int row_index, int new_Ac) {
    // Updates the row in the tables vector at the index w. new_Ac value
    table[row_index].update_Ac(new_Ac);
}

const std::vector<int>& Table::get_rowcol(int row_index) {
    return table[row_index].coloring;
}

int Table::get_row_col_len(int row_index) {
    return table[row_index].coloring.size();
}

int Table::get_vertex_col_index(int v) {
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

void Table::remove_from_rows_coloring(int row_index, int v_index) {
    table[row_index].remove_from_coloring(v_index);
}

int Table::get_rows_key(int row_index) {
    return table[row_index].get_key();
}

int Table::create_row(int coloring) {
    // creates a new row and adds it to the table/table_lookups
    Row r;        
    if(coloring!=-1) r.append_coloring(coloring);
    int index = insert_row(r);
    
    return index;
}

int Table::create_row(int copy_row_index, int coloring) {   
    Row r = table[copy_row_index];
    if(coloring!=-1) r.append_coloring(coloring);
    int index = insert_row(r);

    return index;
}

int Table::copyin_row(Table* child_table, int child_row_index) {
    Row r = child_table->table[child_row_index];
    int index = insert_row(r);

    return index;
}

int Table::copyin_forgetrow(Table* child_table, int child_row_index, int v_index) {    
    Row r(child_table->table[child_row_index]);
    r.remove_from_coloring(v_index);
    r.set_childr_table_ind(child_row_index);    
    int index = insert_row(r);
    return index;
}

int Table::lookup_Ac(int key) {
    /*
     * Given a key, returns the Ac value of the row corresponding
     * to that key.
     */
    int index = lookup_table_index(key);
    return table[index].get_Ac();
}

int Table::get_Ac(int row_index) {
    return table[row_index].get_Ac();
}

int Table::get_rows_childl_table_ind(int row_index) {
    return table[row_index].get_childl_table_ind();
}

int Table::get_rows_childr_table_ind(int row_index) {
    return table[row_index].get_childr_table_ind();
}

void Table::update_rows_childl_table_ind(int row_index, int new_childind) {
    table[row_index].set_childl_table_ind(new_childind);
}

void Table::update_rows_childr_table_ind(int row_index, int new_childind) {
    table[row_index].set_childr_table_ind(new_childind);
}

int Table::insert_row(Row r) {
    // returns the index of the inserted row
    // only inserts if r not already inserted 

    int index = -1;
    if(lookup_table_index(r.get_key()) == -1) {
        index = table.size(); 
        table.push_back(r);       //NOTE makes copy
            
        // key is the concatenated coloring, value is the index in table.
        table_lookups_insert(r.get_key(), index);
    }
    return index;
}

void Table::delete_row(int index) {
    table.erase(table.begin()+index, table.begin()+index+1);  //*3^ni
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

void Table::update_row_add(int row_index, int value) {
    /* must update key in table_lookups whenever we add values
     * to the current coloring.
     * 
     * used for updating a coloring which is already in the table lookups. 
     */
    //1. erase
    int index = lookup_table_index(table[row_index].get_key());
    table_lookups_remove(table[row_index].get_key());
    
    //2. re-insert into table_lookups
    table[row_index].append_coloring(value);
    table_lookups_insert(table[row_index].get_key(), index);
}

int Table::get_table_size() {
    return table.size();
}


