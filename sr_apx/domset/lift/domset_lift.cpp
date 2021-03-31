
#include "sr_apx/domset/lift/domset_lift.hpp"
#include "sr_apx/domset/apx/domset_apx.hpp"

namespace sr_apx {
namespace domset {
namespace lift {

Set naive_lift(const Graph& graph, const Set& editset, const Set& partial) {
	Set domset(editset);
	domset.insert(partial.begin(), partial.end());

	return domset;
}

Set greedy_lift(const Graph& graph, const Set& editset, const Set& partial) {
	Set undominated;
	Map<Set>::const_iterator iu = graph.begin();
	for ( ; iu != graph.end(); ++iu) {
		if (partial.contains(iu->first)) {
			continue;
		}

		undominated.insert(iu->first);
		for (int v : iu->second) {
			if (partial.contains(v)) {
				undominated.erase(iu->first);
				break;
			}
		}
	}

	Set domset = apx::greedy_apx(graph, editset, undominated);
	domset.insert(partial.begin(), partial.end());

	return domset;
}

}}}
