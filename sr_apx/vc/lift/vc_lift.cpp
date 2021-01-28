
#include "sr_apx/vc/lift/vc_lift.hpp"
#include "sr_apx/vc/apx/vc_apx.hpp"
#include "sr_apx/vc/exact/vc_exact.hpp"
#include "sr_apx/bipartite/bipartite.hpp"

namespace sr_apx {
namespace vc {
namespace lift {

Set naive_lift(const Graph& graph, const Set& octset, const Set& partial) {
	Set cover(octset);
	cover.insert(partial.begin(), partial.end());

	return cover;
}

Set greedy_lift(const Graph& graph, const Set& octset, const Set& partial) {
	Set cover(partial);

	Set processed;
	for (Map<Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu) {
		int u = iu->first;
		if (!octset.contains(u)) {
			processed.insert(u);
		}
	}

	for (int u : octset) {
		for (int v : graph.neighbors(u)) {
			if (processed.contains(v) && !cover.contains(v)) {
				cover.insert(u);
				break;
			}
		}

		processed.insert(u);
	}

	return cover;
}

Set apx_lift(const Graph& graph, const Set& octset, const Set& partial) {
	Set subgraph_vertices;
	for (Map<Set>::const_iterator it = graph.begin(); it != graph.end(); ++it) {
		if (octset.contains(it->first) || !partial.contains(it->first)) {
			subgraph_vertices.insert(it->first);
		}
	}

	Graph h = graph.subgraph(subgraph_vertices);

	Set cover = apx::std_apx(h);
	cover.insert(partial.begin(), partial.end());

	return cover;
}

Set oct_lift(const Graph& graph, const Set& octset, const Set& partial) {
	Graph h = graph.subgraph(octset);
	Set cover = apx::std_apx(h);
	cover.insert(partial.begin(), partial.end());

	Set subgraph_vertices;
	for (Map<Set>::const_iterator it = graph.begin(); it != graph.end(); ++it) {
		if (!cover.contains(it->first)) {
			subgraph_vertices.insert(it->first);
		}
	}

	h = graph.subgraph(subgraph_vertices);
	Set bipcover = exact::bip_exact(h);
	cover.insert(bipcover.begin(), bipcover.end());

	return cover;
}

Set bip_lift(const Graph& graph, const Set& octset, const Set& partial) {
	Graph h;
	for (Map<Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu) {
		int u = iu->first;
		if (!octset.contains(u)) {
			continue;
		}

		for (int v : graph.neighbors(u)) {
			if (!octset.contains(v) && !partial.contains(v)) {
				h.add_edge(u, v);
			}
		}
	}

	Set cover = exact::bip_exact(h);

	Set subgraph_vertices;
	for (int u : octset) {
		if (!cover.contains(u)) {
			subgraph_vertices.insert(u);
		}
	}

	h = graph.subgraph(subgraph_vertices);
	Set octcover = apx::std_apx(h);

	cover.insert(partial.begin(), partial.end());
	cover.insert(octcover.begin(), octcover.end());

	return cover;
}

Set recursive_lift(const Graph& graph, const Set& octset, const Set& partial) {
	Set subgraph_vertices;
	for (Map<Set>::const_iterator it = graph.begin(); it != graph.end(); ++it) {
		if (octset.contains(it->first) || !partial.contains(it->first)) {
			subgraph_vertices.insert(it->first);
		}
	}

	Graph h = graph.subgraph(subgraph_vertices);

	Set octset2 = bipartite::vertex_delete(h);
	subgraph_vertices.clear();
	for (Map<Set>::const_iterator it = h.begin(); it != h.end(); ++it) {
		if (!octset2.contains(it->first)) {
			subgraph_vertices.insert(it->first);
		}
	}

	Graph g = h.subgraph(subgraph_vertices);
	Set cover = exact::bip_exact(g);

	cover.insert(partial.begin(), partial.end());
	cover.insert(octset2.begin(), octset2.end());

	return cover;
}

Set recursive_oct_lift(const Graph& graph, const Set& octset, const Set& partial) {
	Graph h = graph.subgraph(octset);
	Set octset2 = bipartite::vertex_delete(h);
	Set subgraph_vertices;
	for (Map<Set>::const_iterator it = h.begin(); it != h.end(); ++it) {
		if (!octset2.contains(it->first)) {
			subgraph_vertices.insert(it->first);
		}
	}

	Graph g = h.subgraph(subgraph_vertices);
	Set cover = exact::bip_exact(g);

	cover.insert(octset2.begin(), octset2.end());
	cover.insert(partial.begin(), partial.end());

	subgraph_vertices.clear();
	for (Map<Set>::const_iterator it = graph.begin(); it != graph.end(); ++it) {
		if (!cover.contains(it->first)) {
			subgraph_vertices.insert(it->first);
		}
	}

	g = graph.subgraph(subgraph_vertices);
	Set bipcover = exact::bip_exact(g);

	cover.insert(bipcover.begin(), bipcover.end());

	return cover;
}

Set recursive_bip_lift(const Graph& graph, const Set& octset, const Set& partial) {
	Graph h;
	for (Map<Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu) {
		int u = iu->first;
		if (!octset.contains(u)) {
			continue;
		}

		for (int v : graph.neighbors(u)) {
			if (!octset.contains(v) && !partial.contains(v)) {
				h.add_edge(u, v);
			}
		}
	}

	Set cover = exact::bip_exact(h);

	Set subgraph_vertices;
	for (int u : octset) {
		if (!cover.contains(u)) {
			subgraph_vertices.insert(u);
		}
	}

	h = graph.subgraph(subgraph_vertices);

	Set octset2 = bipartite::vertex_delete(h);
	subgraph_vertices.clear();
	for (Map<Set>::const_iterator it = h.begin(); it != h.end(); ++it) {
		if (!octset2.contains(it->first)) {
			subgraph_vertices.insert(it->first);
		}
	}

	Graph g = h.subgraph(subgraph_vertices);
	Set octcover = exact::bip_exact(g);

	cover.insert(octcover.begin(), octcover.end());
	cover.insert(octset2.begin(), octset2.end());
	cover.insert(partial.begin(), partial.end());

	return cover;
}

}}}
