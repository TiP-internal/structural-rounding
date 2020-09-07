
#ifndef GRAPH_H
#define GRAPH_H

#include "sr_apx/setmap/setmap.hpp"

namespace sr_apx {

class Graph {
public:
	Map<Set> adjlist;

	Graph() {};
	Graph(int);
	~Graph();
	void add_edge(int, int);
	int size();
	int degree(int);
	bool adjacent(int, int);

	Map<Set>::iterator begin();
	Map<Set>::iterator end();
	Set* neighbors(int);

	Graph* subgraph(Set*);
	void remove_vertex(int);
};

Graph* read_sparse6(const char*);
Graph* read_edge_list(const char*);
Graph* read_dimacs(const char*);

}

#endif
