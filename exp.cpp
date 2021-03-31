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

  printf("graph,read time,k,n,m,den,es,esden,esdenbw,es size,tw (sep),tw (gd),tw dif,sep time,gd time,speedup");
  //printf("graph,read time,n,m,tw (sep),tw (gd),tw dif,sep time,gd time,speedup,");
  printf(",log size,log time");
  //printf("graph name,n,m,read time,tw,log size,log time");
  printf(",edit2 size,edit2 time,edit2 tw (sep),edit2 tw (gd),tw dif,sep time,gd time,speedup,partial2 size (sep),partial2 size (gd)");
  //printf("size dif,sep time,gd time,speedup,lift time (sep),lift time (gd),total2 size (sep),total2 size (gd),size dif,sep time,gd time,speedup");
  printf(",size dif,sep time,gd time,speedup,lift time (sep),lift time (gd),sep time,gd time,speedup");
  //printf(",ed3,ed3 time,ed3 tw,par3,par3 time,tot3,tot3 time");
  //printf(",ed4,ed4 time,ed4 tw,par4,par4 time,tot4,tot4 time");
  //printf(",ed5,ed5 time,ed5 tw,par5,par5 time,tot5,tot5 time");
  printf("\n");

  for (std::vector<std::string>::iterator graph_files_it = graph_files.begin(); graph_files_it != graph_files.end(); graph_files_it++) {
    std::string filename = *graph_files_it;

    /* verify name */
    if (filename.find(".txt") == std::string::npos)
      continue;
    printf("%s,", (filename.substr(0, filename.length()-3)).c_str());

    /* read graph */
    clock_t start = clock();
    sr_apx::Graph graph = sr_apx::read_edge_list((filepath + filename).c_str());
    clock_t end = clock();
    printf("%.2fs,", (double)(end-start)/1000000);

    /* print k (synthetic) */
    int l = filename.find_first_of("k")+1;
    int r = filename.substr(l).find_first_of("_");
    printf("%d,", std::stoi(filename.substr(l, r)));

    /* find n,m */
    int deg = 0;
    for (sr_apx::Map<sr_apx::Set>::const_iterator iu = graph.begin(); iu != graph.end(); ++iu)
      deg += iu->second.size();
    printf("%d,", graph.size());
    printf("%d,", deg >> 1);

    /* synthetic graph properties */
    // print den
    l = filename.find_first_of("r")+1;
    r = filename.substr(l).find_first_of("_");
    printf("%.3f,", std::stod(filename.substr(l, r)));
    // print es
    int pad = l+1;
    l = filename.substr(pad).find_first_of("s")+1;
    r = filename.substr(pad).substr(l).find_first_of("_");
    printf("%d,", std::stoi(filename.substr(pad).substr(l, r)));
    // print esd
    pad += l+r;
    l = filename.substr(pad).find_first_of("d")+1;
    r = filename.substr(pad).substr(l).find_first_of("_");
    printf("%.2f,", std::stod(filename.substr(pad).substr(l, r)));
    // print esdbw
    pad += l+r;
    l = filename.substr(pad).find_first_of("w")+1;
    r = filename.substr(pad).substr(l).find_first_of("_");
    printf("%.2f,", std::stod(filename.substr(pad).substr(l, r)));
    // print es size
    pad += l+r;
    l = filename.substr(pad).find_last_of("e")+1;
    r = filename.substr(pad).substr(l).find_first_of("_");
    printf("%.2f,", std::stod(filename.substr(pad).substr(l, r)));

    /* initial decompositions */
    // separator-based
    start = clock();
    sr_apx::treewidth::Decomposition init_decomp(graph);
    end = clock();
    double sep = (double)(end-start)/1000000;
    int sep_tw = init_decomp.treewidth();
    printf("%d,", sep_tw);
    // greedy degree
    start = clock();
    //sr_apx::treewidth::Decomposition gd_decomp(false);
    //gd_decomp.build_gd_decompisition(graph);
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
    printf("%.3f,", (double)gd_tw/sep_tw);
    // comparisons
    printf("%.2fs,", sep);
    printf("%.2fs,", gd);
    //printf("%.2f," gd/gd_old);
    printf("%.2f,", sep/gd);

    /* log */
    start = clock();
    sr_apx::Set domset = sr_apx::domset::apx::greedy_apx(graph);
    end = clock();
    printf("%d,", domset.size());
    printf("%.2fs", (double)(end-start)/1000000);
    check_ds(graph, domset);


    for (int i = 2; i <= 2; i++) {
      /* edit */
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
      end = clock();
      double edit_time = (double)(end-start)/1000000;
      printf(",%d,%.4f,", edit.size(), edit_time); // edit2 size, edit2 time

      /* edit decomp */
      // separator-based
      start = clock();
      sr_apx::treewidth::Decomposition edit_decomp(sub_g);
      end = clock();
      double sep_edit = (double)(end-start)/1000000;
      int sep_edit_tw = edit_decomp.treewidth();
      printf("%d,", sep_edit_tw); // edit2 tw (sep)
      // greedy degree
      sr_apx::treewidth::Decomposition edit_decomp2(sub_g); // remove later
      start = clock();
      edit_decomp2.build_gd_decomposition(sub_g);
      end = clock();
      double gd_edit = (double)(end-start)/1000000;
      int gd_edit_tw = edit_decomp.treewidth();
      // greedy degree old (for accurate tw - remove after the above line works correctly)
      std::vector<int> gd_edit_order = sr_apx::treewidth::greedy_degree(graph, graph.size());
      sr_apx::Graph gd_edit_plus = sr_apx::treewidth::fill(graph, graph.size(), gd_order);
      edit_decomp2.build_decomposition(gd_edit_plus, gd_edit_order);
      gd_edit_tw = edit_decomp2.treewidth();
      printf("%d,", gd_edit_tw); // edit2 tw (gd)
      // comparisons
      printf("%.3f,", (double)gd_edit_tw/sep_edit_tw); // tw dif (<= better)
      printf("%.2fs,", sep_edit); // sep time
      printf("%.2fs,", gd_edit); // gd time
      printf("%.2f,", sep_edit/gd_edit); // speedup

      /* solve */
      // separator-based
      start = clock();
      sr_apx::Set partial;
      //sr_apx::domset::exact::calculate(sub_g, edit_decomp, partial, opt, sr_apx::domset::exact::Variant::Dom_Set, true);
      end = clock();
      double sep_partial = (double)(end-start)/1000000;
      int sep_partial_sz = partial.size();
      printf("%d,", sep_partial_sz); // partial2 size (sep)
      // greedy degree
      start = clock();
      sr_apx::Set partial2;
      //sr_apx::domset::exact::calculate(sub_g, edit_decomp2, partial2, opt, sr_apx::domset::exact::Variant::Dom_Set, true);
      end = clock();
      double gd_partial = (double)(end-start)/1000000;
      int gd_partial_sz = partial2.size();
      printf("%d,", gd_partial_sz); // partial2 size (gd)
      // comparisons
      printf("%.3f,", (double)gd_partial_sz/sep_partial_sz); // size dif (>= better?)
      printf("%.2fs,", sep_partial); // sep time
      printf("%.2fs,", gd_partial); // gd time
      printf("%.2f,", sep_partial/gd_partial); // speedup

      /* lift */
      // separator-based
      start = clock();
      //sr_apx::Set domset = sr_apx::domset::lift::greedy_lift(graph, edit, partial);
      end = clock();
      double sep_lift = (double)(end-start)/1000000;
      printf("%.2fs,", sep_lift); // lift time (sep)
      // greedy degree
      start = clock();
      //sr_apx::Set domset2 = sr_apx::domset::lift::greedy_lift(graph, edit, partial2);
      end = clock();
      double gd_lift = (double)(end-start)/1000000;
      printf("%.2fs,", gd_lift); // lift time (gd)
      // comparisons
      //printf("%d,", domset.size()); // total2 size (sep)
      //printf("%d,", domset2.size()); // total2 size (gd)
      //printf("%.3f,", (double)domset2.size()/domset.size()); // size dif (>= better?)
      double total_sep = edit_time + sep_edit + sep_partial + sep_lift;
      double total_gd = edit_time + gd_edit + gd_partial + gd_lift;
      printf("%.2fs,", total_sep); // sep time
      printf("%.2fs,", total_gd); // gd time
      printf("%.2f", total_sep/total_gd); // speedup

      //check_ds(graph, domset);
      //check_ds(graph, domset2);
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
