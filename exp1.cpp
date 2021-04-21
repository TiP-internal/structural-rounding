#include <dirent.h>
#include <iostream>
#include <vector>
#include <time.h>

#include "sr_apx/sr_apx.hpp"

void read_directory(const std::string&, std::vector<std::string>&);
void check_ds(const sr_apx::Graph&, const sr_apx::Set&);


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

    printf("graph,read time,tw,n,m,drop,ess,esd,bd,");
    printf("log,log time,");
    //printf("ed8,ed8 time,par8,par8 time,tot8,tot8 time,");
    printf("ed7,ed7 time,par7,par7 time,tot7,tot7 time,");
    printf("ed6,ed6 time,par6,par6 time,tot6,tot6 time,");
    printf("ed5,ed5 time,par5,par5 time,tot5,tot5 time,");
    printf("ed4,ed4 time,par4,par4 time,tot4,tot4 time,");
    printf("ed3,ed3 time,par3,par3 time,tot3,tot3 time,");
    printf("ed2,ed2 time,par2,par2 time,tot2,tot2 time");
    printf("\n");

    for (std::vector<std::string>::iterator graph_files_it = graph_files.begin(); graph_files_it != graph_files.end(); graph_files_it++) {
        std::string filename = *graph_files_it;
    
        /* verify name */
        if (filename.find(".s6") == std::string::npos)
            continue;
        printf("%s,", (filename.substr(0, filename.length()-4)).c_str()); // graph

        /* read graph */
        clock_t start = clock();
        sr_apx::Graph graph = sr_apx::read_sparse6((filepath + filename).c_str());
        clock_t end = clock();
	printf("%.2fs,", (double)(end-start)/1000000); // read time

        /* print tw (synthetic) */
        int l = filename.find_first_of("w")+1;
        int r = filename.substr(l).find_first_of("_");
        int tw = std::stoi(filename.substr(l,r));
        printf("%d,", tw); // tw

        /* print n, m */
        int deg = 0;
        for (sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu)
            deg += iu->second.size();
        printf("%d,", graph.size()); // n
        printf("%d,", deg >> 1); // m

        /* synthetic graph properties */
        // print drop
        l = filename.find_first_of("p")+1;
        r = filename.substr(l).find_first_of("_");
	printf("%.1f,", std::stod(filename.substr(l, r))); // drop
        // print ess
        int pad = l+1;
        l = filename.substr(pad).find_first_of("s")+2;
        r = filename.substr(pad).substr(l).find_first_of("_");
        printf("%.2f,", std::stod(filename.substr(pad).substr(l, r))); // ess
        // print esd
        pad += l+r;
        l = filename.substr(pad).find_first_of("d")+1;
        r = filename.substr(pad).substr(l).find_first_of("_");
	printf("%.3f,", std::stod(filename.substr(pad).substr(l, r))); // esd
        // print bd
        pad += l+r;
        l = filename.substr(pad).find_first_of("d")+1;
        r = filename.substr(pad).substr(l).find_first_of("_");
        printf("%.3f,", std::stod(filename.substr(pad).substr(l, r))); // bd

        /* log */
        start = clock();
        sr_apx::Set domset = sr_apx::domset::apx::greedy_apx(graph);
        end = clock();
        printf("%d,", domset.size()); // log
        printf("%.2fs,", (double)(end-start)/1000000); // log time
        check_ds(graph, domset);

	// fill tw vals we don't edit to
	for (int i = tw; i < 7; i++)
	  printf(",,,,,,");

        for (int i = tw; i > 1; i--) {      
            /* edit */
            start = clock();
	    sr_apx::Set edit = sr_apx::treewidth::vertex_gd_delete(graph, i);
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
            end = clock();
            double edit_time = (double)(end-start)/1000000;
            printf("%d,", edit.size()); // ed
	    printf("%.2fs,", edit_time); // ed time

            /* solve */
	    start = clock();
	    sr_apx::treewidth::Decomposition decomp(true);
	    decomp.build_gd_decomposition(sub_g);
	    //printf("\nbefore solve\n");
	    sr_apx::Set partial = sr_apx::domset::exact::tw_exact(sub_g, decomp, opt);
            end = clock();
	    //printf("after solve\n");
	    printf("%d,", partial.size()); // par
	    double partial_time = (double)(end-start)/1000000;
	    printf("%.2fs,", partial_time); // par time

            /* lift */
            start = clock();
	    sr_apx::Set domset = sr_apx::domset::lift::greedy_lift(graph, edit, partial);
            end = clock();
	    printf("%d,", domset.size()); // tot
            double lift_time = (double)(end-start)/1000000;
	    printf("%.2fs", edit_time + partial_time + lift_time); // tot time
            if (i != 2) printf(",");
	    check_ds(graph, domset);
        }

        printf("\n");
    }

    return 0;
}


void read_directory(const std::string& name, std::vector<std::string>& v) {
    DIR* dirp = opendir(name.c_str());
    struct dirent* dp;
    while ((dp = readdir(dirp)) != NULL)
        v.push_back(dp->d_name);

    closedir(dirp);
}


void check_ds(const sr_apx::Graph& graph, const sr_apx::Set& ds) {
    sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin();
    for ( ; iu != graph.end(); ++iu) {
        if (ds.contains(iu->first))
             continue;

        bool flag = false;
        for (int nbr : iu->second) {
            if (ds.contains(nbr)) {
	        flag = true;
	        break;
            }
        }

        if (!flag)
            printf("%d is not dominated\n", iu->first);
    }
}
