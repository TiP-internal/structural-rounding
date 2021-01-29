
#include "sr_apx/graph/graph.hpp"
#include "sr_apx/util/util.hpp"

#include <stdexcept>
#include <iostream>
#include <fstream>

#include <cstring>  // for strcmp in read_dimacs

namespace sr_apx {

Graph::Graph(int n): adjlist(n) {}

bool Graph::adjacent(int u, int v) const {
	Map<Set>::const_iterator loc = adjlist.find(u);
	if (loc == adjlist.end()) {
		return false;
	}

    return loc->second.contains(v);
}

bool Graph::contains_vertex(int u) const {
	return adjlist.contains(u);
}

void Graph::add_edge(int u, int v) {
	adjlist[u].insert(v);
	adjlist[v].insert(u);
}

void Graph::remove_edge(int u, int v) {
	Map<Set>::iterator loc = adjlist.find(u);
	if (loc == adjlist.end()) {
		return;
	}

	loc->second.erase(v);

	loc = adjlist.find(v);
	if (loc == adjlist.end()) {
		return;
	}

	loc->second.erase(u);
}

int Graph::size() const {
	return adjlist.size();
}

Map<Set>::iterator Graph::begin() {
	return adjlist.begin();
}

Map<Set>::iterator Graph::end() {
	return adjlist.end();
}

Map<Set>::const_iterator Graph::begin() const {
	return adjlist.cbegin();
}

Map<Set>::const_iterator Graph::end()const  {
	return adjlist.cend();
}

Set& Graph::neighbors(int u) {
	return adjlist.at(u);
}

const Set& Graph::neighbors(int u) const {
	return adjlist.at(u);
}

int Graph::degree(int u) const {
	Map<Set>::const_iterator loc = adjlist.find(u);
	if (loc == adjlist.end()) {
		return 0;
	}

	return loc->second.size();
}

Graph Graph::subgraph(const Set& vertices) const {
	Graph subg(vertices.size());

	for (int u : vertices) {
		Map<Set>::const_iterator loc = adjlist.find(u);
		if (loc == adjlist.end()) {
			continue;
		}

		subg.adjlist.insert({u, Set()});

		for (int v : loc->second) {
			if (vertices.contains(v)) {
				subg.add_edge(u, v);
			}
		}
	}

	return subg;
}

void Graph::remove_vertex(int vertex) {
	Map<Set>::iterator loc = adjlist.find(vertex);
	if (loc == adjlist.end()) {
		return;
	}

	for (int nbr : loc->second) {
		adjlist.at(nbr).erase(vertex);
	}

	adjlist.erase(vertex);
}

std::vector<Set> Graph::connected_components() const {
	std::vector<Set> components;

	Set comp;
	std::vector<int> stack;
	Set visited;

	Map<Set>::const_iterator vitr = adjlist.begin();

	while (!stack.empty() || visited.size() < adjlist.size()) {
		int current;
		if (stack.empty()) {
			if (!comp.empty()) {
				components.push_back(std::move(comp));
				comp = Set();
			}

			while (visited.contains(vitr->first)) {
				++vitr;
			}

			current = vitr->first;
			visited.insert(current);
			comp.insert(current);
		}
		else {
			current = stack.back();
			stack.pop_back();
		}

		for (int nbr : adjlist.at(current)) {
			if (!comp.contains(nbr)) {
				stack.push_back(nbr);
				comp.insert(nbr);
				visited.insert(nbr);
			}
		}
	}

	if (!comp.empty()) {
		components.push_back(std::move(comp));
	}

	return components;
}

Graph read_sparse6(const char* filename) {
	std::ifstream f;
	f.open(filename, std::ios::in | std::ios::binary);
	char c[7];
	f.read(c, 1);
	if (c[0] != ':') {
		throw std::invalid_argument("graph is not sparse6");
	}

	int n;
	f.read(c, 1);
	if (c[0] - 63 < 63) {
		n = c[0] - 63;
	}
	else {
		f.read(c, 3);
		if (c[0] - 63 < 63) {
			n = (c[0] - 63) << 12;
			n += (c[1] - 63) << 6;
			n += (c[2] - 63);
		}
		else {
			f.read(&c[3], 4);
			n = (c[1] - 63) << 30;
			n += (c[2] - 63) << 24;
			n += (c[3] - 63) << 18;
			n += (c[4] - 63) << 12;
			n += (c[5] - 63) << 6;
			n += (c[6] - 63);
		}
	}

	static constexpr int buffer_size = 1024;

	int k = sr_apx::util::log2(n);
	Graph graph(n);

	int bitbuffer = 0;
	int bitavailable = 0;

	int position = 0;
	int available = 0;
	char buffer[buffer_size];

	int v = 0;
	while (position < available || !f.eof()) {
		if (bitavailable < 1) {
			if (position == available) {
				f.read(buffer, buffer_size);
				available = f.gcount();
				position = 0;
			}

			bitbuffer = buffer[position++] - 63;
			bitavailable = 6;
		}
		bitavailable -= 1;
		int b = bitbuffer >> bitavailable;
		bitbuffer &= (1 << bitavailable) - 1;

		while (bitavailable < k) {
			if (position == available) {
				f.read(buffer, buffer_size);
				available = f.gcount();
				position = 0;
			}

			bitbuffer = (bitbuffer << 6) + (buffer[position++] - 63);
			bitavailable += 6;
		}

		bitavailable -= k;
		int x = bitbuffer >> bitavailable;
		bitbuffer &= (1 << bitavailable) - 1;

		if (b == 1) {
			v += 1;
		}

		if (x >= n || v >= n) {
			break;
		}
		else if (x > v) {
			v = x;
		}
		else {
			graph.add_edge(x, v);
		}
	}

	f.close();
	return graph;
}

Graph read_edge_list(const char* filename) {
	std::ifstream f;
	f.open(filename, std::ios::in);

	Graph g;

	char s[100];
	f.getline(s, 100);
	while (f.getline(s, 100)) {
		int u, v;
		sscanf(s, "%d %d", &u, &v);
		g.add_edge(u, v);
	}

	f.close();
	return g;
}

Graph read_dimacs(const char* filename) {
    /* Reads graphs in the DIMACS graph format.
     * 
     * Example:
        c comment 1
        c comment 2
        c ...
        p nodes edges 
        e 1 2
        e 2 3
        ....
    */
    
    std::ifstream f;
    f.open(filename, std::ios::in);
    
    Graph g;

    char s[100];
    while (f.getline(s, 100)) {
        
        char type[2];
        int u, v;
        sscanf(s, "%s %d %d", type, &u, &v);
        
        char edge_type[] = "e";
        if (strcmp(type, edge_type) == 0) {
            g.add_edge(u, v);
        } 
    }

    f.close();
    return g;
}

}
