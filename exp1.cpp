#include <iostream>
#include <time.h>

#include "sr_apx/sr_apx.hpp"


void print_graph(const sr_apx::Graph&);

/* pass a partial k-tree edge list file for the argument */
int main(int argc, char* argv[]) {
    std::string file = argv[1];
    sr_apx::Graph g = sr_apx::read_edge_list(file.c_str());
    int n = g.size();

    double total = 0;
    printf("\nfile: %s (n: %d) (k: ?)\n", file.c_str(), n);
    clock_t start = clock();
    std::vector<int> ordering0 = sr_apx::treewidth::greedy_degree(g, n);
    clock_t end = clock();
    total += (double)(end-start);
    printf("greedy_degree  = %.4f%s", (double)(end-start)/1000000, "s\n");

    start = clock();
    sr_apx::Graph g_plus0 = sr_apx::treewidth::fill(g, n, ordering0);
    end = clock();
    total += (double)(end-start);
    printf("fill           = %.4f%s", (double)(end-start)/1000000, "s\n");

    start = clock();
    std::reverse(ordering0.begin(), ordering0.end());
    sr_apx::treewidth::Decomposition d0(g);
    d0.build_decomposition(g, ordering0);
    end = clock();
    total += (double)(end-start);
    printf("tree decomp    = %.4f%s", (double)(end-start)/1000000, "s\n");
    printf("total          = %.4f%s", total/1000000, "s\n\n");



    total = 0;
    start = clock();
    std::vector<int> ordering1 = sr_apx::treewidth::greedy_fill_in(g, n);
    end = clock();
    total += (double)(end-start);
    printf("greedy_fill_in = %.4f%s", (double)(end-start)/1000000, "s\n");

    start = clock();
    sr_apx::Graph g_plus1 = sr_apx::treewidth::fill(g, n, ordering1);
    end = clock();
    total += (double)(end-start);
    printf("fill           = %.4f%s", (double)(end-start)/1000000, "s\n");

    start = clock();
    sr_apx::treewidth::Decomposition d1(g);
    d1.build_decomposition(g, ordering1);
    end = clock();
    total += (double)(end-start);
    printf("tree decomp    = %.4f%s", (double)(end-start)/1000000, "s\n");

    start = clock();
    sr_apx::Graph tm = sr_apx::treewidth::minimal_triangulation(g_plus1);
    end = clock();
    total += (double)(end-start);
    printf("tri minimiz    = %.4f%s", (double)(end-start)/1000000, "s\n");
    printf("total          = %.4f%s", total/1000000, "s\n\n");
}


void print_graph(const sr_apx::Graph& g) {
    sr_apx::Map<sr_apx::Set>::const_iterator iu = g.begin();
    for ( ; iu != g.end(); ++iu)
        printf("v: %d (%d)\n", iu->first, g.degree(iu->first));
}
