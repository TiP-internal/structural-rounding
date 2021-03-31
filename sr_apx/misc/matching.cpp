
#include "sr_apx/misc/matching.hpp"

#include <vector>
#include <deque>

namespace sr_apx {
namespace misc {

Map<int> bipartite_matching(const Graph& graph, const Set& left, const Set& right) {
	Map<int> match;

	bool update = true;

	while (update) {
		update = false;

		std::deque<int> queue;
		Map<int> distance;
		Set unmatched;

		for (int u : left) {
			if (!match.contains(u)) {
				distance[u] = 0;
				queue.push_back(u);
				unmatched.insert(u);
			}
		}

		while (!queue.empty()) {
			int current = queue.front();
			queue.pop_front();

			for (int nbr : graph.neighbors(current)) {
				if (!match.contains(nbr)) {
					continue;
				}

				int m = match.at(nbr);

				if (!distance.contains(m)) {
					distance[m] = distance[current] + 1;
					queue.push_back(m);
				}
			}
		}

		std::vector<int> stack;
		Map<int> parent;

		while (unmatched.size() > 0 || !stack.empty()) {
			if (stack.empty()) {
				stack.push_back(-1);
				int u = *(unmatched.begin());
				unmatched.erase(u);
				stack.push_back(u);
			}

			int current = stack.back();
			stack.pop_back();
			int previous = stack.back();
			stack.pop_back();

			parent[current] = previous;
			for (int nbr : graph.neighbors(current)) {
				if (!match.contains(nbr)) {
					update = true;

					stack.clear();
					int c = current;
					int n = nbr;

					while (c != -1) {
						match[n] = c;
						int temp = -1;
						if (match.contains(c)) {
							temp = match.at(c);
						}
						match[c] = n;
						n = temp;
						c = parent[c];
					}

					break;
				}

				if (distance[match.at(nbr)] == distance[current] + 1) {
					stack.push_back(current);
					stack.push_back(match.at(nbr));
				}
			}

			distance.erase(current);
		}
	}

	return match;
}

}}
