
#include "sr_apx/vc/kernel/lp_kernel.hpp"
#include "sr_apx/vc/exact/vc_exact.hpp"

namespace sr_apx::vc::kernel {

Set** lp_kernel(Graph* g) {

    int n = g->size();
    Graph* h = new Graph();
    Set* left = new Set();
    Set* right = new Set();

    for (auto iu = g->begin(); iu != g->end(); ++iu) {
        int u = iu->first;

        left->insert(u);
        right->insert(u + n);

        for (auto inbr = g->neighbors(u)->begin(); inbr != g->neighbors(u)->end(); ++inbr) {
            int nbr = *inbr;
            h->add_edge(u, nbr + 2 * n);
        }
    }

    Set* cover = exact::bip_exact(h);

    Set* in = new Set();
    Set* out = new Set();

    for (auto iu = g->begin(); iu != g->end(); ++iu) {
        int u = iu->first;
        if (cover->contains(u) && cover->contains(u + 2 * n)) {
            in->insert(u);
        }
        else if (!cover->contains(u) && !cover->contains(u + 2 * n)) {
            out->insert(u);
        }
    }

    delete h;
    delete left;
    delete right;
    delete cover;

    Set** ret = new Set*[2];
    ret[0] = in;
    ret[1] = out;
    return ret;
}

}
