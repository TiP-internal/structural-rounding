
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
public:
    //NOTE could possiby delete the vector after table is finished creating?
    std::vector<int> coloring;      //all possible colorings for the xi verts
    Set* domset_verts = new Set();  //Holds domset soln verts
    int A_c;                        //domset size
    int key;
    std::string key_str = "";
    
    Row() {
        this->A_c = INF;
    }
    
    Row(const Row* r) {  //copy constructor
        for(int i=0; i<r->coloring.size(); i++) {
            this->coloring.push_back(r->coloring[i]);
        }
        this->A_c = r->A_c;
        this->key = r->key;
        this->key_str = r->key_str;

        for(auto it=r->domset_verts->begin(); it!=r->domset_verts->end(); it++) {
            this->domset_verts->insert(*it);
        }
    }
    
    ~Row() {
        delete domset_verts;
    }
    
    void append_coloring(int val) {
        coloring.push_back(val);
        
        //val should only ever be a single int (1,2, or 3)
        std::string k = std::to_string(val); 
        key_str = key_str+k;
        
        //Complexity: "Unspecified, but generally linear in the number 
        //of characters interpreted."
        key = stoi(key_str);   
    }
    
    void remove_from_coloring(int v_index) {
        coloring.erase(coloring.begin()+v_index);  //erase the coloring for v
        
        key_str.erase(v_index, 1); //deletes one characters at index?
        key = stoi(key_str);
    }

};


class Table {
public:
    int label;  //stores the current bag index
    //vector of the vertices in the current bag.
    //table[i].coloring[j] corresponds to the coloring of vertices[j] vertex. 
    //could probably delete after table creation since we have pointer to bag.
    std::vector<int> vertices; 
    
    //Keeping this because we need to iterate over rows in order. 
    std::vector<Row*> table;   //vector containing pointers to row structs
    Map<int> table_lookups;    //key=the unique coloring, value=index of row in the table. 
    
    //Public Functions
    Table(); 
    Table(int);
    ~Table();
    
    void insert_row(Row*);
    void update_row_add(Row*, int);
    int get_vertex_col_index(int);
    int lookup(int);
};


#endif
