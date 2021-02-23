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

    // printf("graph,n,m,tw,read time,decomp time,log size,log time");
    // printf(",edit2 size,edit2 time,edit2 actualtw,partial2 size,partial2 time,total2 size,total2 time");
    // printf(",edit3 size,edit3 time,edit3 actualtw,partial3 size,partial3 time,total3 size,total3 time");
    // printf(",edit4 size,edit4 time,edit4 actualtw,partial4 size,partial4 time,total4 size,total4 time");
    // printf(",edit5 size,edit5 time,edit5 actualtw,partial5 size,partial5 time,total5 size,total5 time");
    // printf("\n");

    printf("graph,read time,k,n,m,den,es,esden,esdenbw,es size,tw (sep),tw (gd),sep time,gd time,speedup,");
    //printf("graph,read time,n,m,tw (sep),tw (gd),sep time,gd time,speedup,");
    //printf("graph,read time,n,m,tw (sep),tw (gd),tw (gfi),sep time,gd time,gfi time,");
    printf("log size,log time\n");

    std::vector<std::string> names;
    std::vector<int> n;
    std::vector<int> m;
    std::vector<double> read_times;
    std::vector<int> tw_sep;
    std::vector<int> tw_gd;
    std::vector<int> tw_gfi;
    std::vector<double> tw_sep_times;
    std::vector<double> tw_gd_times;
    std::vector<double> tw_gfi_times;
    std::vector<int> log_sizes;
    std::vector<double> log_times;
    
    std::vector<int> edits;
    std::vector<double> edit_times;
    std::vector<int> edits_tw_sep;
    // std::vector<std::vector<int>> edits;
    // std::vector<std::vector<double>> edit_times;

    for (std::vector<std::string>::iterator graph_files_it = graph_files.begin(); graph_files_it != graph_files.end(); graph_files_it++) {
        /* verify name */
	std::string filename = *graph_files_it;
        if (filename.find(".txt") == std::string::npos)
            continue;
	names.push_back((filename.substr(0, filename.length()-4)).c_str());
	printf("%s,",(filename.substr(0, filename.length()-4)).c_str());

        /* read graph */
	clock_t start = clock();
	sr_apx::Graph graph = sr_apx::read_edge_list((filepath + filename).c_str());
	clock_t end = clock();
	read_times.push_back((double)(end-start)/1000000);
	printf("%.4fs,", (double)(end-start)/1000000);

	/* print k */
	int l = filename.find_first_of("k")+1;
	int r = filename.substr(l).find_first_of("_");
	printf("%d,", std::stoi(filename.substr(l, r)));

	/* find m */
	int deg = 0;
	for (sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu)
	    deg += iu->second.size();
	n.push_back(graph.size());
	m.push_back(deg >> 1);
	printf("%d,", graph.size());
	printf("%d,", deg >> 1);

	// print den 
	l = filename.find_first_of("r")+1;
	r = filename.substr(l).find_first_of("_");
	printf("%.3f,", std::stod(filename.substr(l, r)));
	/* print es */ 
	int pad = l+1;
	l = filename.substr(pad).find_first_of("s")+1;
	r = filename.substr(pad).substr(l).find_first_of("_");
	printf("%d,", std::stoi(filename.substr(pad).substr(l, r)));
	/* print esd */
	pad += l+r;
	l = filename.substr(pad).find_first_of("d")+1;
	r = filename.substr(pad).substr(l).find_first_of("_");
	printf("%.2f,", std::stod(filename.substr(pad).substr(l, r)));
	/* print esdbw */
	pad += l+r;
	l = filename.substr(pad).find_first_of("w")+1;
	r = filename.substr(pad).substr(l).find_first_of("_");
	printf("%.2f,", std::stod(filename.substr(pad).substr(l, r)));
	/* print es size */
	pad += l+r;
	l = filename.substr(pad).find_last_of("e")+1;
	r = filename.substr(pad).substr(l).find_first_of("_");
	printf("%.2f,", std::stod(filename.substr(pad).substr(l, r)));

	/* initial decompositions */
	// seperator
	start = clock();
	sr_apx::treewidth::Decomposition init_decomp(graph);
	end = clock();
	double sep = (double)(end-start)/1000000;
	tw_sep.push_back(init_decomp.treewidth());
	tw_sep_times.push_back((double)(end-start)/1000000);
	printf("%d,", init_decomp.treewidth());
	// greedy degre
	start = clock();
	sr_apx::treewidth::Decomposition gd_decomp(false);
	gd_decomp.build_gd_decomposition(graph);
	//init_decomp.build_gd_decomposition(graph);
	end = clock();
	double gd = (double)(end-start)/1000000;
	printf("%d,", gd_decomp.treewidth());
	// greedy degree old
	start = clock();
	std::vector<int> gd_order = sr_apx::treewidth::greedy_degree(graph, graph.size());
	sr_apx::Graph gd_plus = sr_apx::treewidth::fill(graph, graph.size(), gd_order);
	init_decomp.build_decomposition(gd_plus, gd_order);
	end = clock();
	double gd_old = (double)(end-start)/1000000;
	tw_gd.push_back(init_decomp.treewidth());
	tw_gd_times.push_back((double)(end-start)/1000000);
	//printf("%d,", init_decomp.treewidth());
	
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
	printf("%.2f,", gd/sep);
	//printf("%.4fs (%.2fx),", gfi, gfi/sep);

	/* log */
	start = clock();
	sr_apx::Set domset = sr_apx::domset::apx::greedy_apx(graph);
	end = clock();
	log_sizes.push_back(domset.size());
	log_times.push_back((double)(end-start)/1000000);
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
	    edits.push_back(edit.size());
	    edit_times.push_back((double)(end-start)/1000000);

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
            // tw_gd.push_back(init_decomp.treewidth());
	    // tw_gd_times.push_back((double)(end-start)/1000000);


            sr_apx::Set domset = sr_apx::domset::lift::greedy_lift(graph, edit, partial);
            printf(",%d,%.4f", domset.size(), time1 + time2);
            check_ds(graph, domset);
	}
    }

    /* print results */
    /* printf("name,");
    for (int i = 0; i < names.size(); i++) {
        printf("%s", names[i].c_str());
        if (i != names.size()-1)
            printf(",");
    }
    printf("\nread time,");
    for (int i = 0; i < read_times.size(); i++) {
        printf("%.4fs", read_times[i]);
        if (i != read_times.size()-1)
            printf(",");
    }
    printf("\nn,");
    for (int i = 0; i < n.size(); i++) {
        printf("%d", n[i]);
        if (i != n.size()-1)
            printf(",");
    }
    printf("\nm,");
    for (int i = 0; i < m.size(); i++) {
	printf("%d", m[i]);
	if (i != m.size()-1)
	    printf(",");
    }
    printf("\ntw (sep),");
    for (int i = 0; i < tw_sep.size(); i++) {
        printf("%d", tw_sep[i]);
        if (i != tw_sep.size()-1)
            printf(",");
    }
    printf("\ntw (gd),");
    for (int i = 0; i < tw_gd.size(); i++) {
        printf("%d", tw_gd[i]);
        if (i != tw_gd.size()-1)
	    printf(",");
    }
    printf("\ntw (gfi),");
    for (int i = 0; i < tw_gfi.size(); i++) {
        printf("%d", tw_gfi[i]);
        if (i != tw_gfi.size()-1)
	    printf(",");
    }
    printf("\nsep time,");
    for (int i = 0; i < tw_sep_times.size(); i++) {
        printf("%.4fs", tw_sep_times[i]);
        if (i != tw_sep_times.size()-1)
            printf(",");
    }
    printf("\ngd time,");
    for (int i = 0; i < tw_gd_times.size(); i++) {
        printf("%.4fs (%.2fx)", tw_gd_times[i], tw_gd_times[i]/tw_sep_times[i]);
        if (i != tw_gd_times.size()-1)
            printf(",");
    }
    printf("\ngfi time,");
    for (int i = 0; i < tw_gfi_times.size(); i++) {
        printf("%.4fs (%.2fx)", tw_gfi_times[i], tw_gfi_times[i]/tw_sep_times[i]);
	if (i != tw_gfi_times.size()-1)
	    printf(",");
    }
    printf("\nlog size,");
    for (int i = 0; i < log_sizes.size(); i++) {
        printf("%d", log_sizes[i]);
        if (i != log_sizes.size()-1)
      	    printf(",");
    }
    printf("\nlog time,");
    for (int i = 0; i < log_times.size(); i++) {
	printf("%.4fs", log_times[i]);
	if (i != log_times.size()-1)
	    printf(",");
    }
    printf("\nedit2 size,");
    for (int i = 0; i < edits.size(); i++) {
        printf("%d", edits[i]);
        if (i != edits.size()-1)
	    printf(",");
    }
    printf("\nedit2 tw,");
    for (int i = 0; i < edits_tw_sep.size(); i++) {
        printf("%d", edits_tw_sep[i]);
        if (i != edits_tw_sep.size()-1)
	    printf(",");
    }
    printf("\nedit2 time,");
    for (int i = 0; i < edit_times.size(); i++) {
        printf("%.4fs", edit_times[i]);
        if (i != edit_times.size()-1)
	    printf(",");
    } */

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
