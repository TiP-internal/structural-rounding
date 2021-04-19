
#include "sr_apx/graph/graph.hpp"
#include "sr_apx/util/util.hpp"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <random>

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

void writebits(long& bitbuffer, int& used, int k, int value, std::ofstream& f) {
    used += k;
    bitbuffer = (bitbuffer << k) + value;

    while (used >= 6) {
        used -= 6;
		f << (char) ((bitbuffer >> used) + 63);
        bitbuffer &= (1 << used) - 1;
    }
}

void write_sparse6(const Graph& graph, const char* filename) {
	std::vector<int> vertices;
	vertices.reserve(graph.size());
	Map<Set>::const_iterator iu = graph.begin();
	for ( ; iu != graph.end(); ++iu) {
		vertices.push_back(iu->first);
	}

	std::sort(vertices.begin(), vertices.end());
	int n = vertices[vertices.size() - 1] + 1;
	int k = util::log2(n);
	int mask = (1 << 6) - 1;

	std::ofstream f(filename, std::ios::out | std::ios::binary);

	f << ':';
	if (n <= 62) {
		f << (char) (63 + n);
	}
	else if (n <= 258047) {
		f << (char) 126;
		f << (char) (((n >> 12) & mask) + 63);
		f << (char) (((n >> 6) & mask) + 63);
		f << (char) ((n & mask) + 63);
	}
	else {
		f << (char) 126 << (char) 126;
		f << (char) (((n >> 30) & mask) + 63);
		f << (char) (((n >> 24) & mask) + 63);
		f << (char) (((n >> 18) & mask) + 63);
		f << (char) (((n >> 12) & mask) + 63);
		f << (char) (((n >> 6) & mask) + 63);
		f << (char) ((n & mask) + 63);
	}

	long bitbuffer = 0;
	int used = 0;
	int prev = 0;

	for (int v : vertices) {
		for (int nbr : graph.neighbors(v)) {
			if (nbr < v) {
				if (v == prev + 1) {
					writebits(bitbuffer, used, 1, 1, f);
				}
				else if (v > prev + 1) {
					writebits(bitbuffer, used, 1, 0, f);
					writebits(bitbuffer, used, k, v, f);
					writebits(bitbuffer, used, 1, 0, f);
				}
				else {
					writebits(bitbuffer, used, 1, 0, f);
				}

				writebits(bitbuffer, used, k, nbr, f);
				prev = v;
			}
		}
	}

	writebits(bitbuffer, used, 1, 1, f);
	writebits(bitbuffer, used, 6 - used, 0, f);

	f << '\n';

	f.close();
}

void write_edge_list(const Graph& graph, const char* filename) {
	std::ofstream f(filename, std::ios::out);

	f << graph.size() << '\n';

	Map<Set>::const_iterator iu = graph.begin();
	for ( ; iu != graph.end(); ++iu) {
		int u = iu->first;
		for (int nbr : iu->second) {
			if (nbr > u) {
				f << u << ' ' << nbr << '\n';
			}
		}
	}

	f.close();
}

Graph shuffle_vertices(const Graph& graph) {
	std::vector<int> vertices;
	vertices.reserve(graph.size());

	Map<Set>::const_iterator iu = graph.begin();
	for (iu = graph.begin(); iu != graph.end(); ++iu) {
		vertices.push_back(iu->first);
	}

	std::mt19937 gen(70801032);
	std::shuffle(vertices.begin(), vertices.end(), gen);

	Map<int> labels;
	labels.reserve(graph.size());

	for (int i = 0; i < vertices.size(); ++i) {
		labels[vertices[i]] = i;
	}

	Graph g;
	for (iu = graph.begin(); iu != graph.end(); ++iu) {
		for (int nbr : iu->second) {
			if (nbr < iu->first) {
				g.add_edge(labels[nbr], labels[iu->first]);
			}
		}
	}

	return g;
}

}
