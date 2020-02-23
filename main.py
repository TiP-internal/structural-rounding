
import cProfile

import sys, os, random
from time import time
from csv import DictWriter

from sr_apx.graph import Graph, read_edge_list, read_sparse6
from sr_apx.setmap import Set
from sr_apx.octset import prescribed_octset, find_octset, verify_bip

from sr_apx.vc.apx import dfs_apx, std_apx, heuristic_apx
from sr_apx.vc.exact import bip_exact
from sr_apx.vc.lift import naive_lift, greedy_lift

def run_apx(apx, graph, n):
    times = []
    sols = []
    for i in range(n):
        random.seed(i)
        start = time()
        cover = apx(graph)
        end = time()
        times.append(end - start)
        sols.append(len(cover))

    avgtime = round(sum(times) / n, 4)
    minsol = min(sols)
    maxsol = max(sols)
    return avgtime, minsol, maxsol

def run_lift(lift, graph, n, octset, partial):
    times = []
    sols = []
    for i in range(n):
        random.seed(i)
        start = time()
        cover = lift(graph, octset, partial)
        end = time()
        times.append(end - start)
        sols.append(len(cover))

    avgtime = round(sum(times) / n, 4)
    minsol = min(sols)
    maxsol = max(sols)
    return avgtime, minsol, maxsol

def main():
    config_args= []
    if (len(sys.argv) > 2):
        config_args = read_config(sys.argv[2])

    filepath = sys.argv[1]
    directory = True;
    if filepath[len(filepath)-1] != '/':
        if (filepath[len(filepath)-3] == '.' and filepath[len(filepath)-2] == 's'
            and filepath[len(filepath)-1] == '6'):
            directory = False;
        else:
            filepath += "/";
    n = 1

    results_dir = os.path.join(os.getcwd(),"results")
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    with open("results/results.csv", "w") as f:
        header = ["name","n","m","dfs time","dfs size","heuristic time","heuristic size","std time","std size","stdrev time","stdrev size","oct size","partial","bip time","naive time","naive size","apx time","apx size","greedy time","greedy size","octfirst time","octfirst size","octfirst break","bipfirst time","bipfirst size","bipfirst break","rec time","rec size","rec break","recoct time","recoct size","recoct break","recbip time","recbip size","recbip break"]
        results = DictWriter(f, header)
        results.writeheader()

        graph_files = []
        if directory:
            if os.access(filepath, os.W_OK):
                for filename in os.listdir(filepath):
                    if filename.endswith(".s6"):
                        graph_files.append(filename)
        else:
            if not "/" in filepath:
                graph_files.append(filepath)
                filepath = ""
            else:
                graph_files.append(filepath[filepath.rfind("/")+1:len(filepath)])
                filepath = filepath[0:filepath.rfind("/")+1]

        for filename in graph_files:
            res = {}

            graphname = filename[0:len(filename)-3]
            print(graphname)
            res["name"] = graphname

            start = time()
            graph = read_sparse6("{}{}".format(filepath, filename))
            end = time()
            print("n: {}".format(len(graph)))
            print("time: {}".format(round(end - start, 4)))
            res["n"] = len(graph)

            t, minsol, maxsol = run_apx(heuristic_apx, graph, n)
            print("heuristic apx")
            print("\tavg time: {}".format(t))
            print("\tmin size: {}".format(minsol))
            print("\tmax size: {}".format(maxsol))
            res["heuristic time"] = t
            res["heuristic size"] = minsol

            t, minsol, maxsol = run_apx(dfs_apx, graph, n)
            print("dfs apx")
            print("\tavg time: {}".format(t))
            print("\tmin size: {}".format(minsol))
            print("\tmax size: {}".format(maxsol))
            res["dfs time"] = t
            res["dfs size"] = minsol

            t, minsol, maxsol = run_apx(std_apx, graph, n)
            print("std apx")
            print("\tavg time: {}".format(t))
            print("\tmin size: {}".format(minsol))
            res["std time"] = t
            res["std size"] = minsol

            # left, right, octset = prescibed_octset(graph, "{}{}.oct".format(filepath, graphname))

            start = time()
            left, right, octset = find_octset(graph)
            # left, right, octset = verify_bip(graph, set())

            bippart = Set()
            for v in left:
                bippart.add(v)
            for v in right:
                bippart.add(v)


            partial = bip_exact(graph.subgraph(bippart))
            end = time()

            print("bip solve")
            print("\tavg time: {}".format(round(end - start, 4)))

            print(len(partial))

            res["oct size"] = len(octset)
            res["partial"] = len(partial)
            res["bip time"] = round(end - start, 4)

            t, minsol, maxsol = run_lift(naive_lift, graph, n, octset, partial)
            print("naive lift")
            print("\tavg time: {}".format(t))
            print("\tmin size: {}".format(minsol))
            print("\tmax size: {}".format(maxsol))
            res["naive time"] = t
            res["naive size"] = minsol

            t, minsol, maxsol = run_lift(greedy_lift, graph, n, octset, partial)
            print("greedy lift")
            print("\tavg time: {}".format(t))
            print("\tmin size: {}".format(minsol))
            print("\tmax size: {}".format(maxsol))
            res["greedy time"] = t
            res["greedy size"] = minsol

            results.writerow(res)
            del graph

            print()

def read_config(name):
    if not os.path.isfile(name):
        print("config file error")
        exit(1)

    args = []
    with open(name) as infile:
        for line in infile:
            if len(line) > 0:
                args.append(line[0:line.rfind(" ")])

    infile.close()
    return args

if __name__ == "__main__":
    # cProfile.run("main()")
    main()
