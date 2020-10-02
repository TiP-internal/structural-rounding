
#include <dirent.h>
//#include <fstream>
#include <iostream>
#include <vector>

#include "time.h"

#include "graph.hpp"
#include "bipartite.hpp"
#include "setmap.hpp"

#include "vc_apx.hpp"
#include "vc_exact.hpp"
#include "vc_lift.hpp"

#include "treewidth.hpp"
#include "domset_exact.hpp"
#include "domset_apx.hpp"

// helper function declarations ////////////////////////////////////////////////

double sum(double* vals, int len);
int min(int* vals, int len);
int max(int* vals, int len);
void read_directory(const std::string& name, std::vector<std::string>& v);

// main.py functions ///////////////////////////////////////////////////////////

double run_apx(Set* (*apx)(Graph*), Graph* graph, int n, int &minsol, int &maxsol) {
	double* times = new double[n];
	int* sols = new int[n];
	for (int i = 0; i < n; i++) {
		clock_t start = clock();
		Set* cover = apx(graph);
		clock_t end = clock();
		times[i] = (double) (end - start);
		sols[i] = cover->size();
		delete cover;
	}

	double avgtime = sum(times, n) / n;
	minsol = min(sols, n);
	maxsol = max(sols, n);
	delete[] times;
	delete[] sols;
	return avgtime;
}

double run_lift(Set* (*lift)(Graph*, Set*, Set*), Graph* graph, int n, Set* octset,
				Set* partial, int &minsol, int &maxsol) {
	double* times = new double[n];
	int* sols = new int[n];
	for (int i = 0; i < n; i++) {
		clock_t start = clock();
		Set* cover = lift(graph, octset, partial);
		clock_t end = clock();
		times[i] = (double) end - start;
		sols[i] = cover->size();
		delete cover;
	}

	double avgtime = sum(times, n) / n;
	minsol = min(sols, n);
	maxsol = max(sols, n);
	delete[] times;
	delete[] sols;
	return avgtime;
}

int main(int argc, char* argv[]) {
	std::string filepath = argv[1];
	bool directory = true;
	if (filepath[filepath.size()-1] != '/') {
		if (filepath[filepath.size()-3] == '.' && filepath[filepath.size()-2] == 's'
		    && filepath[filepath.size()-1] == '6')
			directory = false;
		else
			filepath += "/";
	}
	int n = 1;

	std::vector<std::string> graph_files;
	if (directory)
    	read_directory(filepath, graph_files);
	else {
		if (filepath.find("/") == std::string::npos) {
			graph_files.push_back(filepath);
			filepath = "";
		}
		else {
			graph_files.push_back(filepath.substr(filepath.rfind("/")+1, filepath.length()-2));
			filepath = filepath.substr(0, filepath.rfind("/")+1);
		}
	}

	printf("graph name,n,m,read time,tw,log size,log time,");
	printf("edit2 size,edit2 time,edit2 actualtw,partial2 size,partial2 time,total2 size,total2 time,");
	printf("edit3 size,edit3 time,edit3 actualtw,partial3 size,partial3 time,total3 size,total3 time,");
	printf("edit4 size,edit4 time,edit4 actualtw,partial4 size,partial4 time,total4 size,total4 time,");
	printf("edit5 size,edit5 time,edit5 actualtw,partial5 size,partial5 time,total5 size,total5 time\n");

    for (std::vector<std::string>::iterator graph_files_it = graph_files.begin();
		graph_files_it != graph_files.end(); graph_files_it++) {
		std::string filename = *graph_files_it;
		if (filename.find(".s6") == std::string::npos)
			continue;

		printf("%s,", (filename.substr(0, filename.length()-3)).c_str());

		clock_t start = clock();
		Graph* graph = read_sparse6((filepath + filename).c_str());
		clock_t end = clock();
		printf("%d,", graph->size());
		int deg = 0;
		for (int u : *graph) {
			deg += graph->degree(u);
		}
		printf("%d,", deg >> 1);
		printf("%.4f,", (double)(end-start)/1000000);

		TreeDecomp init(graph);
		printf("%d,", init.treewidth());

		start = clock();
		Set* domset = logn_apx(graph);
		end = clock();
		printf("%d,", domset->size());
		printf("%.4f", (double)(end-start)/1000000);
		delete domset;

		for (int i = 2; i <= 5; i++) {
			start = clock();
			Set* opt = new Set();
			Set* edit = treewidth_nodeedit(graph, opt, i, true);
			Set* verts = graph->get_vertices();
			Set* rest = verts->set_minus(edit);
			Graph* sub_g = graph->subgraph_wsingles(rest);
			end = clock();
			double time1 = (double)(end-start)/1000000;
			printf(",%d,%.4f,", edit->size(), time1);
			delete verts;
			delete rest;

			start = clock();
			TreeDecomp* decomp = new TreeDecomp(sub_g);
			int domset_size = calc_min_domset(sub_g, decomp, opt, Variant::Dom_Set);
			end = clock();
			double time2 = (double)(end-start)/1000000;
			printf("%d,%d,%.4f,%d,%.4f", decomp->treewidth(), domset_size, time2, domset_size + edit->size(), time1 + time2);
			delete opt;
			delete edit;
			delete sub_g;
			delete decomp;
		}

		printf("\n");

		delete graph;
		graph = NULL;
	}

	return 0;
}

// original helper functions ///////////////////////////////////////////////////

double sum(double* vals, int len) {
	double sum = 0.0;
	for (int i = 0; i < len; i++)
		sum += vals[i];
	return sum;
}

int min(int* vals, int len) {
	int min = vals[0];
	for (int i = 1; i < len; i++)
		if (vals[i] < min)
			min = vals[i];
	return min;
}

int max(int* vals, int len) {
	int max = vals[0];
	for (int i = 1; i < len; i++)
		if (vals[i] < max)
			max = vals[i];
	return max;
}

// borrowed helper functions ///////////////////////////////////////////////////

/* http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html */
void read_directory(const std::string& name, std::vector<std::string>& v) {
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}
