
#include "sr_apx/bipartite/bipartite.hpp"

#include <vector>
#include <deque>
#include <iostream>
#include <fstream>

namespace sr_apx::bipartite {

Set prescribed_octset(const Graph& graph, const char* filename) {
	std::ifstream f;
	f.open(filename, std::ios::in);

	Set octset;

	char s[100];
	while (f.getline(s, 100)) {
		int u;
		sscanf(s, "%d", &u);
		octset.insert(u);
	}

	f.close();
	return octset;
}

std::tuple<Set, Set, Set> verify_bipartite(const Graph& graph, const Set& os) {
	Set left;
	Set right;
	Set octset;

	Set visited(os);

	Map<Set>::const_iterator git = graph.begin();

	std::deque<int> queue;
	int current;

	while (visited.size() < graph.size() || !queue.empty()) {
		if (queue.empty()) {
			while (visited.contains(git->first)) {
				++git;
			}

			current = git->first;
			visited.insert(current);
		}
		else {
			current = queue.front();
			queue.pop_front();
		}

		for (int nbr : graph.neighbors(current)) {
			if (os.contains(nbr) || octset.contains(nbr)) {
				continue;
			}

			if (!visited.contains(nbr)) {
				visited.insert(nbr);
				queue.push_back(nbr);
				continue;
			}

			if (left.contains(nbr)) {
				right.insert(current);
			}
			else if (right.contains(nbr)) {
				left.insert(current);
			}
		}

		if (left.contains(current) && right.contains(current)) {
			left.erase(current);
			right.erase(current);
			octset.insert(current);
		}
		else if (!left.contains(current) && !right.contains(current)) {
			left.insert(current);
		}
	}

	return {std::move(octset), std::move(left), std::move(right)};
}

void remove_indset(const Graph& graph, Set& available) {
	Map<int> deg;
	std::vector<Set> revdeg;

	int maxdeg = 0;
	for (int u : available) {
		int degree = 0;
		for (int v : graph.neighbors(u)) {
			if (available.contains(v)) {
				degree += 1;
			}
		}

		maxdeg = degree > maxdeg ? degree : maxdeg;
		deg[u] = degree;
		revdeg.resize(maxdeg + 1);
		revdeg[degree].insert(u);
	}

	while (deg.size() > 0) {
		int mindeg;
		for (mindeg = 0; revdeg[mindeg].size() == 0; ++mindeg);

		int u = *(revdeg[mindeg].begin());
		revdeg[mindeg].erase(u);
		deg.erase(u);
		available.erase(u);

		for (int v : graph.neighbors(u)) {
			if (!deg.contains(v)) {
				continue;
			}

			revdeg[deg[v]].erase(v);
			deg.erase(v);

			for (int w : graph.neighbors(v)) {
				if (!deg.contains(w)) {
					continue;
				}

				int degree = deg[w];
				revdeg[degree].erase(w);
				deg[w] = degree - 1;
				revdeg[degree - 1].insert(w);
			}
		}
	}
}

Set vertex_delete(const Graph& graph) {
	Set octset;
	for (Map<Set>::const_iterator it = graph.begin(); it != graph.end(); ++it) {
		octset.insert(it->first);
	}

	remove_indset(graph, octset);
	remove_indset(graph, octset);

	return octset;
}

} // end of sr_apx namepsace
