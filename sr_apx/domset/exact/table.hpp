
#ifndef TABLE_H
#define TABLE_H

#include <vector> 
#include <string>

#include <limits>  //inf
#include <map>

#include "graph.hpp"
#include "setmap.hpp"
#include "treedecomp.hpp"  //for po_bag struct

using namespace std;


#define IN_DOMSET 1
#define DOMINATED 2
#define NOT_DOMINATED 3

#define INF (int)std::numeric_limits<double>::infinity()

class Row {
    /*
     * A_c: 4 bytes
     * key: 4 bytes
     * coloring: k*4 bytes
     */
private:
    int A_c;                        //domset size
    int key;   
    
    int childl_ind;
    int childr_ind;
public:
    std::vector<int> coloring;      //all possible colorings for the xi verts
    
    Row();
    Row(const Row*);
    ~Row(); 
    
    void append_coloring(int val);
    void remove_from_coloring(int v_index);
    
    int get_Ac();
    void update_Ac(int);
    
    int get_key();
    void update_key(int);
    
    int get_childl_table_ind();
    void set_childl_table_ind(int);
    
    int get_childr_table_ind();
    void set_childr_table_ind(int);
};


class Table {
    /*
     * std::vector<Row*> table:  3^k*4bytes 
     *      - Keeping this (over a Map<Row*>) because we need 
     *        to add rows to the end of the table. ie. in introduce table update. 
     *      - Updating rows as we loop over the table, and appending new rows to
     *        the end. If it was a Map, we could possibly loop over rows we just
     *        added, or miss rows we need to update. 
     * 
     * Map<int> table_lookups: 3^k*2*4bytes 
     *      - This is a map w. format: 
     *          key=the unique coloring, value=index of row in the table.
     *      - Use this since the colorings are unique per table. This is used 
     *        for constant time lookups. 
     *      - Requires a little extra work when updating tables 
     *          ie. we must update the keys for introduce and forget tables. 
     *          This means we must remove the {key, val} from the Map, update key,
     *          then readd to the Map. --Though each step should be constant time. 
     *          (except when updating the forget keys, this takes linear).
     * 
     * std::vector<int> vertices: k*4bytes
     *      - Vector of the vertices in the current bag.
     *      - Each row coloring corresponds to a vertex. Indices match between 
     *        colorings and the vertices. 
     * 
     * (3^k+3^k*2+k)*4bytes
     * 
     */
private:
    std::vector<Row*> table;  
    Map<int> table_lookups;    
    
public:
    std::vector<int> vertices; 
    
    //Public Functions
    Table();
    Table(int); 
    
    ~Table();
    
    void update_row_add(Row*, int);
    void insert_row(Row*);
    void delete_row(int);
    
    void table_lookups_insert(int, int);
    void table_lookups_remove(int); 
    bool table_lookups_contains(int);
    
    int lookup_table_index(int);
    int lookup_Ac(int);
    
    Row* lookup_row(int);
    Row* get_row(int);
    
    Row* create_row(Row*, int);
    Row* create_row(int);
    
    int get_table_size();
    int get_vertex_col_index(int);
    void print_tablelookups();
};


#endif
