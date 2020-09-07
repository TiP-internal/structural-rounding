
#include "sr_apx/vc/lift/vc_lift.hpp"
#include "sr_apx/vc/apx/vc_apx.hpp"
#include "sr_apx/vc/exact/vc_exact.hpp"
#include "sr_apx/bipartite/bipartite.hpp"

namespace sr_apx::vc::lift {

Set* naive_lift(Graph* graph, Set* octset, Set* partial) {
	Set* cover = new Set();
	for (Set::iterator iu = octset->begin(); iu != octset->end(); ++iu) {
		int u = *iu;
		cover->insert(u);
	}

	for (Set::iterator iu = partial->begin(); iu != partial->end(); ++iu) {
		int u = *iu;
		cover->insert(u);
	}

	return cover;
}

Set* greedy_lift(Graph* graph, Set* octset, Set* partial) {
	Set* cover = new Set();
	for (Set::iterator iu = partial->begin(); iu != partial->end(); ++iu) {
		int u = *iu;
		cover->insert(u);
	}

	Set processed;
	for (auto iu = graph->begin(); iu != graph->end(); ++iu) {
		int u = iu->first;
		if (!octset->contains(u)) {
			processed.insert(u);
		}
	}

	for (Set::iterator iu = octset->begin(); iu != octset->end(); ++iu) {
		int u = *iu;
		for (Set::iterator iv = graph->neighbors(u)->begin(); iv != graph->neighbors(u)->end(); ++iv) {
			int v = *iv;
			if (processed.contains(v) && !cover->contains(v)) {
				cover->insert(u);
				break;
			}
		}

		processed.insert(u);
	}

	return cover;
}

Set* apx_lift(Graph* graph, Set* octset, Set* partial) {
	Set* subgraph_vertices = new Set();
	for (auto it = graph->begin(); it != graph->end(); ++it) {
		if (octset->contains(it->first) || !partial->contains(it->first)) {
			subgraph_vertices->insert(it->first);
		}
	}

	Graph* h = graph->subgraph(subgraph_vertices);
	delete subgraph_vertices;

	Set* cover = apx::std_apx(h);
	for (Set::iterator it = partial->begin(); it != partial->end(); ++it) {
		cover->insert(*it);
	}

	delete h;
	return cover;
}

Set* oct_lift(Graph* graph, Set* octset, Set* partial) {
	Graph* h = graph->subgraph(octset);
	Set* cover = apx::std_apx(h);
	for (Set::iterator it = partial->begin(); it != partial->end(); ++it) {
		cover->insert(*it);
	}

	delete h;

	Set* subgraph_vertices = new Set();
	for (auto it = graph->begin(); it != graph->end(); ++it) {
		if (!cover->contains(it->first)) {
			subgraph_vertices->insert(it->first);
		}
	}

	h = graph->subgraph(subgraph_vertices);
	Set* bipcover = exact::bip_exact(h);
	for (Set::iterator it = bipcover->begin(); it != bipcover->end(); ++it) {
		cover->insert(*it);
	}

	delete h;
	delete subgraph_vertices;
	delete bipcover;

	return cover;
}

Set* bip_lift(Graph* graph, Set* octset, Set* partial) {
	Graph* h = new Graph();
	for (auto iu = graph->begin(); iu != graph->end(); ++iu) {
		int u = iu->first;
		if (!octset->contains(u)) {
			continue;
		}

		Set* nbrs = graph->neighbors(u);
		for (Set::iterator iv = nbrs->begin(); iv != nbrs->end(); ++iv) {
			int v = *iv;
			if (!octset->contains(v) && !partial->contains(v)) {
				h->add_edge(u, v);
			}
		}
	}

	Set* cover = exact::bip_exact(h);
	delete h;

	Set* subgraph_vertices = new Set();
	for (Set::iterator it = octset->begin(); it != octset->end(); ++it) {
		if (!cover->contains(*it)) {
			subgraph_vertices->insert(*it);
		}
	}

	h = graph->subgraph(subgraph_vertices);
	delete subgraph_vertices;
	Set* octcover = apx::std_apx(h);
	delete h;

	for (Set::iterator it = partial->begin(); it != partial->end(); ++it) {
		cover->insert(*it);
	}

	for (Set::iterator it = octcover->begin(); it != octcover->end(); ++it) {
		cover->insert(*it);
	}

	delete octcover;

	return cover;
}

Set* recursive_lift(Graph* graph, Set* octset, Set* partial) {
	Set* subgraph_vertices = new Set();
	for (auto it = graph->begin(); it != graph->end(); ++it) {
		if (octset->contains(it->first) || !partial->contains(it->first)) {
			subgraph_vertices->insert(it->first);
		}
	}

	Graph* h = graph->subgraph(subgraph_vertices);
	delete subgraph_vertices;

	Set* octset2 = bipartite::vertex_delete(h);
	subgraph_vertices = new Set();
	for (auto it = h->begin(); it != h->end(); ++it) {
		if (!octset2->contains(it->first)) {
			subgraph_vertices->insert(it->first);
		}
	}

	Graph* g = h->subgraph(subgraph_vertices);
	delete h;
	delete subgraph_vertices;
	Set* cover = exact::bip_exact(g);

	for (Set::iterator it = partial->begin(); it != partial->end(); ++it) {
		cover->insert(*it);
	}

	for (Set::iterator it = octset2->begin(); it != octset2->end(); ++it) {
		cover->insert(*it);
	}

	delete octset2;
	delete g;

	return cover;
}

Set* recursive_oct_lift(Graph* graph, Set* octset, Set* partial) {
	Graph* h = graph->subgraph(octset);
	Set* octset2 = bipartite::vertex_delete(h);
	Set* subgraph_vertices = new Set();
	for (auto it = h->begin(); it != h->end(); ++it) {
		if (!octset2->contains(it->first)) {
			subgraph_vertices->insert(it->first);
		}
	}

	Graph* g = h->subgraph(subgraph_vertices);
	delete h;
	delete subgraph_vertices;
	Set* cover = exact::bip_exact(g);
	delete g;

	for (Set::iterator it = octset2->begin(); it != octset2->end(); ++it) {
		cover->insert(*it);
	}

	delete octset2;

	for (Set::iterator it = partial->begin(); it != partial->end(); ++it) {
		cover->insert(*it);
	}

	subgraph_vertices = new Set();
	for (auto it = graph->begin(); it != graph->end(); ++it) {
		if (!cover->contains(it->first)) {
			subgraph_vertices->insert(it->first);
		}
	}

	g = graph->subgraph(subgraph_vertices);
	delete subgraph_vertices;
	Set* bipcover = exact::bip_exact(g);
	delete g;

	for (Set::iterator it = bipcover->begin(); it != bipcover->end(); ++it) {
		cover->insert(*it);
	}

	delete bipcover;

	return cover;
}

Set* recursive_bip_lift(Graph* graph, Set* octset, Set* partial) {
	Graph* h = new Graph();
	for (auto iu = graph->begin(); iu != graph->end(); ++iu) {
		int u = iu->first;
		if (!octset->contains(u)) {
			continue;
		}

		Set* nbrs = graph->neighbors(u);
		for (Set::iterator iv = nbrs->begin(); iv != nbrs->end(); ++iv) {
			int v = *iv;
			if (!octset->contains(v) && !partial->contains(v)) {
				h->add_edge(u, v);
			}
		}
	}

	Set* cover = exact::bip_exact(h);
	delete h;

	Set* subgraph_vertices = new Set();
	for (Set::iterator it = octset->begin(); it != octset->end(); ++it) {
		if (!cover->contains(*it)) {
			subgraph_vertices->insert(*it);
		}
	}

	h = graph->subgraph(subgraph_vertices);
	delete subgraph_vertices;

	Set* octset2 = bipartite::vertex_delete(h);
	subgraph_vertices = new Set();
	for (auto it = h->begin(); it != h->end(); ++it) {
		if (!octset2->contains(it->first)) {
			subgraph_vertices->insert(it->first);
		}
	}

	Graph* g = h->subgraph(subgraph_vertices);
	delete subgraph_vertices;
	delete h;
	Set* octcover = exact::bip_exact(g);
	delete g;

	for (Set::iterator it = octcover->begin(); it != octcover->end(); ++it) {
		cover->insert(*it);
	}

	delete octcover;

	for (Set::iterator it = octset2->begin(); it != octset2->end(); ++it) {
		cover->insert(*it);
	}

	delete octset2;

	for (Set::iterator it = partial->begin(); it != partial->end(); ++it) {
		cover->insert(*it);
	}

	return cover;
}

}
