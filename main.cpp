
#include <dirent.h>
//#include <fstream>
#include <iostream>
#include <vector>

#include "time.h"

#include "sr_apx/sr_apx.hpp"

// helper function declarations ////////////////////////////////////////////////

double sum(double* vals, int len);
int min(int* vals, int len);
int max(int* vals, int len);
void read_directory(const std::string& name, std::vector<std::string>& v);

// main.py functions ///////////////////////////////////////////////////////////

double run_apx(sr_apx::Set (*apx)(const sr_apx::Graph&), const sr_apx::Graph& graph, int n, int &minsol, int &maxsol) {
	double times[n];
	int sols[n];
	for (int i = 0; i < n; i++) {
		clock_t start = clock();
		sr_apx::Set cover = apx(graph);
		clock_t end = clock();
		times[i] = (double) (end - start);
		sols[i] = cover.size();
	}

	double avgtime = sum(times, n) / n;
	minsol = min(sols, n);
	maxsol = max(sols, n);
	return avgtime;
}

double run_lift(sr_apx::Set (*lift)(const sr_apx::Graph&, const sr_apx::Set&, const sr_apx::Set&), const sr_apx::Graph& graph, int n, const sr_apx::Set& octset, const sr_apx::Set& partial, int &minsol, int &maxsol) {
	double times[n];
	int sols[n];
	for (int i = 0; i < n; i++) {
		clock_t start = clock();
		sr_apx::Set cover = lift(graph, octset, partial);
		clock_t end = clock();
		times[i] = (double) end - start;
		sols[i] = cover.size();
	}

	double avgtime = sum(times, n) / n;
	minsol = min(sols, n);
	maxsol = max(sols, n);
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

    for (std::vector<std::string>::iterator graph_files_it = graph_files.begin();
		graph_files_it != graph_files.end(); graph_files_it++) {
		std::string filename = *graph_files_it;
		if (filename.find(".s6") == std::string::npos)
			continue;

		printf("%s\n", (filename.substr(0, filename.length()-3)).c_str());

		clock_t start = clock();
		sr_apx::Graph graph = sr_apx::read_sparse6((filepath + filename).c_str());
		clock_t end = clock();
		printf("n: %d\n", graph.size());
		printf("time: %.4f\n", (double)(end-start)/1000000);

		int minsol;
		int maxsol;
		double t;

		t = run_apx(sr_apx::vc::apx::heuristic_apx, graph, n, minsol, maxsol);
		printf("heuristic apx\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		t = run_apx(sr_apx::vc::apx::dfs_apx, graph, n, minsol, maxsol);
		printf("dfs apx\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		t = run_apx(sr_apx::vc::apx::std_apx, graph, n, minsol, maxsol);
		printf("std apx\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);

		start = clock();
		sr_apx::Set oct = sr_apx::bipartite::vertex_delete(graph);
		sr_apx::Set left;
		sr_apx::Set right;
		std::tie(std::ignore, left, right) = sr_apx::bipartite::verify_bipartite(graph, oct);

		left.insert(right.begin(), right.end());

		sr_apx::Set partial = sr_apx::vc::exact::bip_exact(graph.subgraph(left));
		end = clock();

		printf("bip solve\n");
		printf("\tavg time: %.4f\n", (double)(end-start)/1000000);

		printf("%d\n", partial.size());

		t = run_lift(sr_apx::vc::lift::naive_lift, graph, n, oct, partial, minsol, maxsol);
		printf("naive lift\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		t = run_lift(sr_apx::vc::lift::greedy_lift, graph, n, oct, partial, minsol, maxsol);
		printf("greedy lift\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		t = run_lift(sr_apx::vc::lift::apx_lift, graph, n, oct, partial, minsol, maxsol);
		printf("apx lift\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		t = run_lift(sr_apx::vc::lift::oct_lift, graph, n, oct, partial, minsol, maxsol);
		printf("oct lift\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		t = run_lift(sr_apx::vc::lift::bip_lift, graph, n, oct, partial, minsol, maxsol);
		printf("bip lift\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		t = run_lift(sr_apx::vc::lift::recursive_lift, graph, n, oct, partial, minsol, maxsol);
		printf("recursive lift\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		t = run_lift(sr_apx::vc::lift::recursive_oct_lift, graph, n, oct, partial, minsol, maxsol);
		printf("recursive oct lift\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		t = run_lift(sr_apx::vc::lift::recursive_bip_lift, graph, n, oct, partial, minsol, maxsol);
		printf("recursive bip lift\n");
		printf("\tavg time: %.4f\n", t/1000000);
		printf("\tmin size: %d\n", minsol);
		printf("\tmax size: %d\n", maxsol);

		printf("\n");
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
