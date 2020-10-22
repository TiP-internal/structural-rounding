
#include "sr_apx/vc/kernel/lp_kernel.hpp"
#include "sr_apx/vc/exact/vc_exact.hpp"

namespace sr_apx::vc::kernel {

std::tuple<Set, Set> lp_kernel(const Graph& g) {
    Graph h(g.size() * 2);
    Set left;
    Set right;

    for (Map<Set>::const_iterator iu = g.begin(); iu != g.end(); ++iu) {
        int u = iu->first;

        left.insert(u);
        right.insert(-u - 1);

        for (int nbr : g.neighbors(u)) {
            h.add_edge(u, -nbr - 1);
        }
    }

    Set cover = exact::bip_exact(h);

    Set in;
    Set out;

    for (Map<Set>::const_iterator iu = g.begin(); iu != g.end(); ++iu) {
        int u = iu->first;
        if (cover.contains(u) && cover.contains(-u - 1)) {
            in.insert(u);
        }
        else if (!cover.contains(u) && !cover.contains(-u - 1)) {
            out.insert(u);
        }
    }

    return {std::move(in), std::move(out)};
}

}
