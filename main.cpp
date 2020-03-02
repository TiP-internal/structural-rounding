
#include <dirent.h>
#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <iterator>

#include "time.h"

#include "graph.hpp"
#include "octset.hpp"
#include "setmap.hpp"

#include "vc_apx.hpp"
#include "vc_exact.hpp"
#include "vc_lift.hpp"


int min(int* vals, int len);
int max(int* vals, int len);
double sum(double* vals, int len);
bool in_array(std::string arg, std::vector<std::string> vec);
std::map<std::string, std::string> parse_config(int argc, char* argv[]);
void read_directory(const std::string& name, std::vector<std::string>& v);


double run_apx(Set* (*apx)(Graph*), Graph* graph, int n, int &minsol, int &maxsol) {
	double* times = new double[n];
	int* sols = new int[n];
	for (int i = 0; i < n; i++) {
		clock_t start = clock();
		Set* cover = apx(graph);
		clock_t end = clock();
		times[i] = (double) (end - start);
		sols[i] = cover->size();
	}

	double avgtime = sum(times, n) / n;
	minsol = min(sols, n);
	maxsol = max(sols, n);
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
	}

	double avgtime = sum(times, n) / n;
	minsol = min(sols, n);
	maxsol = max(sols, n);
	return avgtime;
}


bool is_s6(std::string filepath) {
	if (filepath[filepath.size()-3] == '.' && filepath[filepath.size()-2] == 's'
		&& filepath[filepath.size()-1] == '6')
		return true;
	return false;
}


int main(int argc, char* argv[]) {
	/* parse */
	std::map<std::string, std::string> config_args = parse_config(argc, argv);

	/* format filepath */
	std::string filepath = config_args["graph"];
	bool directory = true;
	if (filepath[filepath.size()-1] != '/') {
		if (is_s6(filepath))
			directory = false;
		else
			filepath += "/";
	}

	/* append graph files */
	std::vector<std::string> graph_files;
	if (directory)
    	read_directory(filepath, graph_files);
	else {
		if (filepath.find("/") == std::string::npos) {
			graph_files.push_back(filepath);
			filepath = "";
		}
		else {
			graph_files.push_back(filepath.substr(filepath.rfind("/")+1, filepath.length()));
			filepath = filepath.substr(0, filepath.rfind("/")+1);
		}
	}

	/* iterate over graph files */
    for (std::vector<std::string>::iterator graph_files_it = graph_files.begin();
		graph_files_it != graph_files.end(); graph_files_it++) {
		std::string filename = *graph_files_it;

		std::string graphname = filename.substr(0, filename.length()-3).c_str();
		printf("%s\n", graphname.c_str());

		clock_t start = clock();
		Graph* graph;
		if (is_s6((filepath + filename).c_str()))
			graph = read_sparse6((filepath + filename).c_str());
		clock_t end = clock();
		printf("n: %d\n", graph->size());
		printf("time: %.4f\n", (double)(end-start)/1000000);

		/* problem */
		if (config_args["problem"] == "vertex_cover") {
			int n = 1;
			int minsol;
			int maxsol;
			double t;

			/* approx sol to problem */
			if (config_args["approx"] != "") {
				if (config_args["approx"] == "heuristic") {
					t = run_apx(heuristic_apx, graph, n, minsol, maxsol);
					printf("heuristic apx\n");
					printf("\tavg time: %.4f\n", t/1000000);
					printf("\tmin size: %d\n", minsol);
					printf("\tmax size: %d\n", maxsol);
				}
				else if (config_args["approx"] == "dfs") {
					t = run_apx(dfs_apx, graph, n, minsol, maxsol);
					printf("dfs apx\n");
					printf("\tavg time: %.4f\n", t/1000000);
					printf("\tmin size: %d\n", minsol);
					printf("\tmax size: %d\n", maxsol);
				}
				else if (config_args["approx"] == "std") {
					t = run_apx(std_apx, graph, n, minsol, maxsol);
					printf("std apx\n");
					printf("\tavg time: %.4f\n", t/1000000);
					printf("\tmin size: %d\n", minsol);
				}
			}

			/* class */
			if (config_args["class"] == "bipartite") {
				/* edit algorithm */
				if (config_args["edit"] == "remove_octset") {
					start = clock();
					OctDecomp* oct = find_octset(graph);

					Set* bippart = new Set();
					for (Set::Iterator left_it = oct->left->begin(); left_it != oct->left->end(); left_it++) {
						int left = *left_it;
						bippart->insert(left);
					}
					for (Set::Iterator right_it = oct->right->begin(); right_it != oct->right->end(); right_it++) {
						int right = *right_it;
						bippart->insert(right);
					}

					Set* partial = bip_exact(graph->subgraph(bippart));
					end = clock();

					printf("bip solve\n");
					printf("\tavg time: %.4f\n", (double)(end-start)/1000000);
					printf("%d\n", partial->size());

					/* lift algorithm */
					if (config_args["lift"] != "") {
						if (config_args["lift"] == "greedy") {
							t = run_lift(naive_lift, graph, n, oct->octset, partial, minsol, maxsol);
							printf("naive lift\n");
							printf("\tavg time: %.4f\n", t/1000000);
							printf("\tmin size: %d\n", minsol);
							printf("\tmax size: %d\n", maxsol);
						}
						else if (config_args["lift"] == "naive") {
							t = run_lift(greedy_lift, graph, n, oct->octset, partial, minsol, maxsol);
							printf("greedy lift\n");
							printf("\tavg time: %.4f\n", t/1000000);
							printf("\tmin size: %d\n", minsol);
							printf("\tmax size: %d\n", maxsol);
						}
					}
				}
			}
			graph = NULL;
			printf("\n");
		}
	}

	return 0;
}


void required_arg_error(std::string argument) {
	std::string arg_flags;
	if (argument == "problem")
		arg_flags = "-p/--" + argument;
	else if (argument == "graph")
		arg_flags = "-g/--" + argument;
	printf("%s%s%s\n" ,"error: argument ", arg_flags.c_str(), ": expected one argument");
	exit(1);
}


void unknown_choice_error(std::string argument, std::string choice, std::vector<std::string> list) {
	std::string arg_flags;
	if (argument == "problem")
		arg_flags = "-p/--" + argument;
	else if (argument == "class")
		arg_flags = "-c/--" + argument;
	else if (argument == "approx")
		arg_flags = "-a/--" + argument;
	else if (argument == "edit")
		arg_flags = "-e/--" + argument;
	else if (argument == "lift")
		arg_flags = "-l/--" + argument;
		
	std::string choices;
	int count = 0;
	for (std::vector<std::string>::iterator it = list.begin(); it != list.end(); it++) {
		std::string choice = *it;
		choices += "\'" + choice + "\'";
		if (count != list.size()-1)
			choices += ", ";
		count ++;
	}

	printf("%s%s%s%s%s%s%s\n", "error: argument ", arg_flags.c_str(), ": invalid choice: \'", choice.c_str(), "\' (choose from ", choices.c_str(), ")");
	exit(1);
}


std::map<std::string, std::string> parse_config(int argc, char* argv[]) {
	std::vector<std::string> problems {"vertex_cover"};
	std::vector<std::string> classes {"bipartite"};
	std::vector<std::string> edits {"remove_octset"};
	std::vector<std::string> lifts {"greedy", "naive"};
	std::vector<std::string> approx {"dfs", "heuristic", "std"};

	std::map<std::string, std::string> config_args;

	while(1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"problem", required_argument, 0, 'p'},
            {"class",   required_argument, 0, 'c'},
            {"edit",    required_argument, 0, 'e'},
            {"lift",    required_argument, 0, 'l'},
            {"approx",  required_argument, 0, 'a'},
            {"graph",   required_argument, 0, 'g'},
			{"spec",    required_argument, 0, 's'},
			{"results", required_argument, 0, 'r'},
            {0,         0,                 0, 0}
		};

		int opt = getopt_long(argc, argv, ":p:c:e:l:a:g:s:r:", long_options, &option_index);
		if (opt == -1)
			break;
		switch(opt) {
			case 'p':
				config_args["problem"] = optarg;
                break;
            case 'c':
				config_args["class"] = optarg;
				break;
			case 'e':
				config_args["edit"] = optarg;
				break;
			case 'l':
				config_args["lift"] = optarg;
				break;
			case 'a':
				config_args["approx"] = optarg;
				break;
            case 'g':
				config_args["graph"] = optarg;
                break;
			case 's':
				config_args["spec"] = optarg;
                break;
			case 'r':
				config_args["results"] = optarg;
                break;
            // case ':':
            //     printf("option needs a value\n");
            //     break;
            // case '?':
            //     printf("unknown option: %c\n", optopt);
            //     break;
        }
	}

	/* check that required arguments exist */
	if (config_args["problem"] == "")
		required_arg_error("problem");
	if (config_args["graph"] == "")
		required_arg_error("graph");
	// (not yet required)
	// if (config_args["results"] == "")
	// 	required_arg_error("results")

	/* check that optional arguments are within their choices */
	if (!(in_array(config_args["problem"], problems)))
		unknown_choice_error("problem", config_args["problem"], problems);
	if (!(in_array(config_args["class"], classes)))
		unknown_choice_error("class", config_args["class"], classes);
	if (!(in_array(config_args["edit"], edits)))
		unknown_choice_error("edit", config_args["edit"], edits);
	if (!(in_array(config_args["lift"], lifts)))
		unknown_choice_error("lift", config_args["lift"], lifts);
	if (!(in_array(config_args["approx"], approx)))
		unknown_choice_error("approx", config_args["approx"], approx);

	return config_args;
}


bool in_array(std::string arg, std::vector<std::string> vec) {
	std::vector<std::string>::iterator it = std::find(vec.begin(), vec.end(), arg);
	std::string match = *it;
	if (arg != match)
		return false;
	return true;
}


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


void read_directory(const std::string& name, std::vector<std::string>& v) {
    DIR* dirp = opendir(name.c_str());
	if (dirp == NULL)
		return;
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
		if (is_s6(dp->d_name))
        	v.push_back(dp->d_name);
    }
    closedir(dirp);
}
