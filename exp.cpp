
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
		if (filename.find(".s6") == std::string::npos)
			continue;

		printf("%s", (filename.substr(0, filename.length()-3)).c_str());

		clock_t start = clock();
		sr_apx::Graph graph = sr_apx::read_sparse6((filepath + filename).c_str());
		clock_t end = clock();
		printf(",%d", graph.size());
		int deg = 0;
		sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin();
		for ( ; iu != graph.end(); ++iu) {
			deg += iu->second.size();
		}
		printf(",%d", deg >> 1);
		printf(",%.4f", (double)(end-start)/1000000);

		sr_apx::treewidth::Decomposition init(graph);
		printf(",%d", init.treewidth());

		// start = clock();
		// Set* domset = logn_apx(graph);
		// end = clock();
		// printf(",%d", domset->size());
		// printf(",%.4f", (double)(end-start)/1000000);
		// delete domset;

		for (int i = 2; i <= 5; i++) {
			start = clock();
			sr_apx::Set edit = sr_apx::treewidth::vertex_delete(graph, i);
			end = clock();
			double time1 = (double)(end-start)/1000000;
			printf(",%d,%.4f", edit.size(), time1);

			// start = clock();
			// TreeDecomp* decomp = new TreeDecomp(sub_g);
			// int domset_size = calculate(sub_g, decomp, NULL, opt, Variant::Dom_Set, false);
			// end = clock();
			// double time2 = (double)(end-start)/1000000;
			// printf(",%d,%d,%.4f,%d,%.4f", decomp->treewidth(), domset_size, time2, domset_size + edit->size(), time1 + time2);
			// delete opt;
			// delete edit;
			// delete sub_g;
			// delete decomp;
		}

		printf("\n");
	}

	return 0;
}
