
#ifndef TABLE_H
#define TABLE_H

#include <vector>
#include <deque>
#include <string>

#include <limits>       // inf
#include <map>
#include <iostream>     // throw
#include <math.h>       // pow

#include "graph.hpp"
#include "setmap.hpp"
#include "treewidth.hpp"

using namespace std;


#define IN_DOMSET 1
#define DOMINATED 2
#define NOT_DOMINATED 3

#define INF (int)std::numeric_limits<double>::infinity()

class Row {
    /*
     * The child table ind functions are what track how the solutions are created.
     * ie are used during the second pass of the algorithm to construct the
     * solution set.
     */
private:
    int A_c;                        //domset size
    unsigned long long int key;     

    int childl_ind;
    int childr_ind;
public:
    std::vector<int> coloring;      //all possible colorings for the xi verts

    Row();
    Row(const Row&);
    Row& operator=(const Row&);
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
     * std::vector<Row*> table:  3^k*bytes
     *      - Keeping this (over a Map<Row*>) because we need
     *        to add rows to the end of the table. ie. in introduce table update.
     *      - Updating rows as we loop over the table, and appending new rows to
     *        the end. If it was a Map, we could possibly loop over rows we just
     *        added, or miss rows we need to update.
     *
     * Map<int> table_lookups: 3^k*2*bytes
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
     * std::vector<int> vertices: k*bytes
     *      - Vector of the vertices in the current bag.
     *      - Each row coloring corresponds to a vertex. Indices match between
     *        colorings and the vertices.
     * 
     * po_bag:
     *      - associated with the current table. used for correct postorder traversal.
     *
     * (3^k+3^k*2+k)*bytes
     *
     */
private:
    std::deque<Row> table;  
    Map<int> table_lookups;    
    po_bag  tables_pobag; 

public:
    std::vector<int> vertices;

    Table();
    Table& operator=(const Table&);
    ~Table();
    
    // Table specific functions
    void set_pobag(po_bag);
    po_bag get_pobag();
    
    int get_table_size();
    
    int create_row(int, int);
    int create_row(int);
    int insert_row(Row);
    int copyin_row(Table*, int);
    int copyin_forgetrow(Table*, int, int);
    int copyin_introrow(Table*, int, int);
    void update_row_add(int, int);
    
    void delete_row(int);
    void pop_front_row();
    
    int get_rows_Ac(int); 
    int get_rows_childl_table_ind(int);
    int get_rows_childr_table_ind(int);
    
    void update_rows_childl_table_ind(int, int);
    void update_rows_childr_table_ind(int, int);
    
    void remove_from_rows_coloring(int, int);
    
    
    // Functions for specific rows in the table vector. 
    void update_Ac(int, int);
    int get_rows_key(int);
    int get_vertex_col_index(int);
    const std::vector<int>& get_rowcol(int);

    void table_lookups_insert(unsigned long long int, int);
    void table_lookups_remove(unsigned long long int);
    bool table_lookups_contains(unsigned long long int);

    int lookup_table_index(unsigned long long int);
    int lookup_Ac(unsigned long long int);
    
    void update_table_lookups(int, int);
    

    //NOTE for testing
    void print_tablelookups() {
        printf("\n===============Table Lookups==================\n");
        for(auto it=table_lookups.begin(); it!=table_lookups.end(); it++) {
            int key = *it;
            printf("coloring=%d, \t table_index=%d\n", key, table_lookups[key]);
        }
        
        printf("\n===========================================\n\n\n");
    }
    
};


#endif
