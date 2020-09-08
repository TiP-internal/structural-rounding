
#include "sr_apx/vc/apx/vc_apx.hpp"

#include <vector>

namespace sr_apx::vc::apx {

Set dfs_apx(const Graph& g) {
	Set cover;

	std::vector<int> stack;
	Set visited;

	Map<Set>::const_iterator vitr = g.begin();

	while (visited.size() < g.size()) {
		int current;
		if (stack.empty()) {
			while (visited.contains(vitr->first)) {
				++vitr;
			}

			current = vitr->first;
		}
		else {
			current = stack.back();
			stack.pop_back();
			int previous = stack.back();
			stack.pop_back();

			if (visited.contains(current)) {
				continue;
			}

			cover.insert(previous);
		}

		visited.insert(current);

		for (int nbr : g.neighbors(current)) {
			if (!visited.contains(nbr)) {
				stack.push_back(current);
				stack.push_back(nbr);
			}
		}
	}


	return cover;
}

void remove_vertex(const Graph& g, Map<int>& deg, Map<Set>& revdeg, int u) {
	for (int nbr : g.neighbors(u)) {
		if (!deg.contains(nbr)) {
			continue;
		}

		int degree = deg[nbr];
		revdeg[degree].erase(nbr);

		if (degree == 1) {
			deg.erase(nbr);
		}
		else {
			deg[nbr] = degree - 1;
			revdeg[degree - 1].insert(nbr);
		}
	}

	deg.erase(u);
}

Set heuristic_apx(const Graph& g) {
	Set cover;

	Map<int> deg(g.size());
	Map<Set> revdeg;
	int maxdeg = 0;

	for (Map<Set>::const_iterator it = g.begin(); it != g.end(); ++it) {
		int degree = g.degree(it->first);
		deg[it->first] = degree;
		revdeg[degree].insert(it->first);
		maxdeg = degree > maxdeg ? degree : maxdeg;
	}

	while (deg.size() > 0) {
		while (revdeg[maxdeg].empty()) {
			--maxdeg;
		}

		int u = *(revdeg[maxdeg].begin());
		revdeg[maxdeg].erase(u);
		cover.insert(u);
		remove_vertex(g, deg, revdeg, u);
	}

	return cover;
}

Set std_apx(const Graph& g) {
	Set cover;

	Map<int> deg(g.size());
	Map<Set> revdeg;
	int maxdeg = 0;

	for (Map<Set>::const_iterator it = g.begin(); it != g.end(); ++it) {
		int degree = g.degree(it->first);
		deg[it->first] = degree;
		revdeg[degree].insert(it->first);
		maxdeg = degree > maxdeg ? degree : maxdeg;
	}

	while (deg.size() > 0) {
		while (revdeg[maxdeg].empty()) {
			--maxdeg;
		}

		int u = *(revdeg[maxdeg].begin());
		int v;
		int md = 0;
		for (int nbr : g.neighbors(u)) {
			if (!deg.contains(nbr)) {
				continue;
			}

			if (g.degree(nbr) > md) {
				md = g.degree(nbr);
				v = nbr;
			}
		}

		revdeg[maxdeg].erase(u);
		cover.insert(u);
		remove_vertex(g, deg, revdeg, u);

		revdeg[deg[v]].erase(v);
		cover.insert(v);
		remove_vertex(g, deg, revdeg, v);
	}

	return cover;
}

}
