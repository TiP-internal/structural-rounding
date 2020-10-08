
#ifndef GRAPH_H
#define GRAPH_H

#include "sr_apx/setmap/setmap.hpp"

#include <vector>

namespace sr_apx {

class Graph {
	Map<Set> adjlist;

public:
	Graph() = default;
	explicit Graph(int);

	void add_edge(int, int);
	void remove_vertex(int);

	int size() const;
	int degree(int) const;
	bool adjacent(int, int) const;
	bool contains_vertex(int) const;

	Map<Set>::iterator begin();
	Map<Set>::iterator end();
	Map<Set>::const_iterator begin() const;
	Map<Set>::const_iterator end() const;

	Set& neighbors(int);
	const Set& neighbors(int) const;

	Graph subgraph(const Set&) const;

	std::vector<Set> connected_components() const;
};

Graph read_sparse6(const char*);
Graph read_edge_list(const char*);
Graph read_dimacs(const char*);

}

#endif
