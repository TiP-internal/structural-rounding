#include <dirent.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <chrono>

#include "sr_apx/sr_apx.hpp"


void read_directory(const std::string&, std::vector<std::string>&);
void check_ds(const sr_apx::Graph&, const sr_apx::Set&);

int main(int argc, char* argv[]) {
    std::string filepath = argv[1];
    bool directory = true;
    if (filepath[filepath.size()-1] != '/') {
        if (filepath[filepath.size()-4] == '.' && filepath[filepath.size()-3] == 't'
	    && filepath[filepath.size()-2] == 'x' && filepath[filepath.size()-1] == 't')
	    directory = false;
        else
	    filepath += "/";
    }

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

    //printf("graph,read time,k,n,m,den,es,esden,esdenbw,es size,tw (sep),tw (gd),tw dif,sep time,gd time,speedup,");
    printf("graph,read time,n,m,tw (sep),tw (gd),tw dif,sep time,gd time,speedup,");
    printf("log size,log time\n");

    for (std::vector<std::string>::iterator graph_files_it = graph_files.begin(); graph_files_it != graph_files.end(); graph_files_it++) {
        /* verify name */
	std::string filename = *graph_files_it;
        if (filename.find(".txt") == std::string::npos)
            continue;
	printf("%s,",(filename.substr(0, filename.length()-4)).c_str());

        /* read graph */
	clock_t start = clock();
	sr_apx::Graph graph = sr_apx::read_edge_list((filepath + filename).c_str());
	clock_t end = clock();
	printf("%.4fs,", (double)(end-start)/1000000);

	/* print k */
	/*int l = filename.find_first_of("k")+1;
	int r = filename.substr(l).find_first_of("_");
	printf("%d,", std::stoi(filename.substr(l, r)));

	/* find m */
	int deg = 0;
	for (sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu)
	    deg += iu->second.size();
	printf("%d,", graph.size());
	printf("%d,", deg >> 1);

	/* print den */ 
	/*l = filename.find_first_of("r")+1;
	r = filename.substr(l).find_first_of("_");
	printf("%.3f,", std::stod(filename.substr(l, r)));
	/* print es */ 
	/*int pad = l+1;
	l = filename.substr(pad).find_first_of("s")+1;
	r = filename.substr(pad).substr(l).find_first_of("_");
	printf("%d,", std::stoi(filename.substr(pad).substr(l, r)));
	/* print esd */
	/*pad += l+r;
	l = filename.substr(pad).find_first_of("d")+1;
	r = filename.substr(pad).substr(l).find_first_of("_");
	printf("%.2f,", std::stod(filename.substr(pad).substr(l, r)));
	/* print esdbw */
	/*pad += l+r;
	l = filename.substr(pad).find_first_of("w")+1;
	r = filename.substr(pad).substr(l).find_first_of("_");
	printf("%.2f,", std::stod(filename.substr(pad).substr(l, r)));
	/* print es size */
	/*pad += l+r;
	l = filename.substr(pad).find_last_of("e")+1;
	r = filename.substr(pad).substr(l).find_first_of("_");
	printf("%.2f,", std::stod(filename.substr(pad).substr(l, r)));

	/* initial decompositions */
	// seperator
	start = clock();
	sr_apx::treewidth::Decomposition init_decomp(graph);
	end = clock();
	double sep = (double)(end-start)/1000000;
	int sep_tw = init_decomp.treewidth();
	printf("%d,", sep_tw);
	
	// greedy degre
	start = clock();
	//sr_apx::treewidth::Decomposition gd_decomp(false);
	//gd_decomp.build_gd_decomposition(graph);
	init_decomp.build_gd_decomposition(graph);
	end = clock();
	double gd = (double)(end-start)/1000000;
	int gd_tw = init_decomp.treewidth();
	//printf("%d,", gd_tw);

	// greedy degree old
	start = clock();
	std::vector<int> gd_order = sr_apx::treewidth::greedy_degree(graph, graph.size());
	sr_apx::Graph gd_plus = sr_apx::treewidth::fill(graph, graph.size(), gd_order);
	init_decomp.build_decomposition(gd_plus, gd_order);
	end = clock();
	double gd_old = (double)(end-start)/1000000;
	gd_tw = init_decomp.treewidth();
	printf("%d,", gd_tw);
	printf("%.3f,", (double) gd_tw/sep_tw);
	
	// greedy fill in
	/*start = clock();
	std::vector<int> gfi_order = sr_apx::treewidth::greedy_fill_in(graph, graph.size());
	sr_apx::Graph gfi_plus = sr_apx::treewidth::fill(graph, graph.size(), gfi_order);
	init_decomp.build_decomposition(gfi_plus, gfi_order);
	end = clock();
	double gfi = (double)(end-start)/1000000;
	tw_gfi.push_back(init_decomp.treewidth());
	tw_gfi_times.push_back((double)(end-start)/1000000);
	printf("%d,", init_decomp.treewidth());*/
	printf("%.4fs,", sep);
	printf("%.4fs,", gd);
	//printf("%.2f,", gd/gd_old);
	printf("%.2f,", sep/gd);
	//printf("%.4fs (%.2fx),", gfi, gfi/sep);

	/* log */
	start = clock();
	sr_apx::Set domset = sr_apx::domset::apx::greedy_apx(graph);
	end = clock();
	printf("%d,", domset.size());
	printf("%.4fs\n",(double)(end-start)/1000000);
	check_ds(graph, domset);

	/* edit */
	for (int i = 2; i <= 1; i++) {
	    start = clock();
	    sr_apx::Set edit = sr_apx::treewidth::vertex_delete(graph, i);
	    sr_apx::Graph sub_g(graph.size() - edit.size());
            sr_apx::Set opt;
            for (sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu) {
                if (edit.contains(iu->first))
                    continue;

                for (int v : iu->second) {
                    if (!edit.contains(v))
                        sub_g.add_edge(iu->first, v);
                }
            }

            for (sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu) {
                if (!edit.contains(iu->first))
                    continue;

                for (int v : iu->second) {
                    if (!edit.contains(v))
                       opt.insert(v);
                }
            }

	    // record edits
	    end = clock();
	    double time1 = (double)(end-start)/1000000;

	    // 
	    start = clock();
            sr_apx::treewidth::Decomposition decomp(sub_g);
            sr_apx::Set partial;
	    sr_apx::domset::exact::calculate(sub_g, decomp, partial, opt, sr_apx::domset::exact::Variant::Dom_Set, true);
            end = clock();
	    double time2 = (double)(end-start)/1000000;
	    printf(",%d,%d,%.4f", decomp.treewidth(), partial.size(), time2);
	    // greedy degre
	    // start = clock();
	    // std::vector<int> gd_order = sr_apx::treewidth::greedy_degree(graph, graph.size());
	    // sr_apx::Graph gd_plus = sr_apx::treewidth::fill(graph, graph.size(), gd_order);
	    // init_decomp.build_decomposition(gd_plus, gd_order);
	    // end = clock();


            sr_apx::Set domset = sr_apx::domset::lift::greedy_lift(graph, edit, partial);
            printf(",%d,%.4f", domset.size(), time1 + time2);
            check_ds(graph, domset);
	}
    }

    return 0;
}

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
