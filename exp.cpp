
#include <dirent.h>
#include <iostream>
#include <vector>

#include "time.h"

#include "sr_apx/sr_apx.hpp"

/* http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html */
void read_directory(const std::string& name, std::vector<std::string>& v) {
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

void check_ds(const sr_apx::Graph& graph, const sr_apx::Set& ds) {
    sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin();
    for ( ; iu != graph.end(); ++iu) {
        if (ds.contains(iu->first)) {
            continue;
        }

        bool flag = false;
        for (int nbr : iu->second) {
            if (ds.contains(nbr)) {
                flag = true;
                break;
            }
        }

        if (!flag) {
            printf("%d is not dominated\n", iu->first);
        }
    }
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

	printf("graph name,n,m,read time,tw,log size,log time");
	printf(",edit2 size,edit2 time,edit2 actualtw,partial2 size,partial2 time,total2 size,total2 time");
	printf(",edit3 size,edit3 time,edit3 actualtw,partial3 size,partial3 time,total3 size,total3 time");
	printf(",edit4 size,edit4 time,edit4 actualtw,partial4 size,partial4 time,total4 size,total4 time");
	printf(",edit5 size,edit5 time,edit5 actualtw,partial5 size,partial5 time,total5 size,total5 time");
	printf("\n");

    for (std::vector<std::string>::iterator graph_files_it = graph_files.begin();
		graph_files_it != graph_files.end(); graph_files_it++) {
		std::string filename = *graph_files_it;
		if (filename.find(".txt") == std::string::npos)
			continue;

		printf("%s", (filename.substr(0, filename.length()-3)).c_str());

		clock_t start = clock();
		sr_apx::Graph graph = sr_apx::read_edge_list((filepath + filename).c_str());
		clock_t end = clock();
		printf(",%d", graph.size());
		int deg = 0;
		sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin();
		for ( ; iu != graph.end(); ++iu) {
			deg += iu->second.size();
		}
		printf(",%d", deg >> 1);
		printf(",%.4f", (double)(end-start)/1000000);

		sr_apx::treewidth::Decomposition init(false);
        init.build_decomposition(graph);
		printf(",%d", init.treewidth());

		start = clock();
		sr_apx::Set domset = sr_apx::domset::apx::greedy_apx(graph);
		end = clock();
		printf(",%d", domset.size());
		printf(",%.4f", (double)(end-start)/1000000);
        check_ds(graph, domset);

		for (int i = 2; i <= 5; i++) {
			start = clock();
			sr_apx::Set edit = sr_apx::treewidth::vertex_delete(graph, i);
            sr_apx::Graph sub_g(graph.size() - edit.size());
            sr_apx::Set opt;
            for (sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu) {
                if (edit.contains(iu->first)) {
                    continue;
                }

                for (int v : iu->second) {
                    if (!edit.contains(v)) {
                        sub_g.add_edge(iu->first, v);
                    }
                }
            }

            for (sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu) {
                if (!edit.contains(iu->first)) {
                    continue;
                }

                for (int v : iu->second) {
                    if (!edit.contains(v)) {
                        opt.insert(v);
                    }
                }
            }

			end = clock();
			double time1 = (double)(end-start)/1000000;
			printf(",%d,%.4f", edit.size(), time1);

			start = clock();
            sr_apx::treewidth::Decomposition decomp(sub_g);
            sr_apx::Set partial = sr_apx::domset::exact::tw_exact(sub_g, decomp, opt);
			end = clock();
			double time2 = (double)(end-start)/1000000;
			printf(",%d,%d,%.4f", decomp.treewidth(), partial.size(), time2);

            sr_apx::Set domset = sr_apx::domset::lift::greedy_lift(graph, edit, partial);

            printf(",%d,%.4f", domset.size(), time1 + time2);
            check_ds(graph, domset);
		}

		printf("\n");
	}

	return 0;
}
