
#include "sr_apx/vc/exact/vc_exact.hpp"
#include "sr_apx/bipartite/bipartite.hpp"
#include "sr_apx/misc/matching.hpp"

#include <stdexcept>
#include <vector>

namespace sr_apx {
namespace vc {
namespace exact {

Set bip_exact(const Graph& graph) {
	Set empty;
	Set octset;
	Set left;
	Set right;
	std::tie(octset, left, right) = bipartite::verify_bipartite(graph, empty);

	if (octset.size() > 0) {
		throw std::invalid_argument("not a bipartite graph");
	}

	Map<int> match = misc::bipartite_matching(graph, left, right);

	Set cover;
	Set visited;
	std::vector<int> stack;

	for (int u : left) {
		cover.insert(u);
		if (!match.contains(u)) {
			visited.insert(u);
			stack.push_back(u);
		}
	}

	int current;
	while (!stack.empty()) {
		current = stack.back();
		stack.pop_back();

		cover.erase(current);

		for (int nbr : graph.neighbors(current)) {
			if (!match.contains(nbr)) {
				cover.insert(nbr);
				continue;
			}

			int m = match.at(nbr);

			if (m != current) {
				cover.insert(nbr);
			}

			if (!visited.contains(m)) {
				visited.insert(m);
				stack.push_back(m);
			}
		}
	}

	return cover;
}

}}}
