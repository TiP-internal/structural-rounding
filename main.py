
import cProfile

import sys, os, argparse, yaml
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
    for _ in range(n):
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
    for _ in range(n):
        start = time()
        cover = lift(graph, octset, partial)
        end = time()
        times.append(end - start)
        sols.append(len(cover))

    avgtime = round(sum(times) / n, 4)
    minsol = min(sols)
    maxsol = max(sols)
    return avgtime, minsol, maxsol


def is_yaml(file):
    return file.lower().endswith('.yaml')


def is_s6(file):
    return file.lower().endswith('.s6')


def main():
    # parse config
    config_args = parse_config()

    # format filepath
    filepath = config_args.graph
    directory = True
    if not filepath.endswith('/'):
        if os.path.isdir(os.path.join(os.getcwd(), filepath)):
            filepath += '/'
        elif os.path.isfile(os.path.join(os.getcwd(), filepath)):
            if is_s6(filepath):
                directory = False
            else:
                print('error: ' + filepath + ' isn\'t a supported graph file type')
                exit(1)
        else:
            print('error: ' + filepath + ' isn\'t a graph file nor a directory')
            exit(1)

    #write to results file
    with open(config_args.results, "w") as f:
        header = ["name","n","m","dfs time","dfs size","heuristic time","heuristic size","std time","std size","stdrev time","stdrev size","oct size","partial","bip time","naive time","naive size","apx time","apx size","greedy time","greedy size","octfirst time","octfirst size","octfirst break","bipfirst time","bipfirst size","bipfirst break","rec time","rec size","rec break","recoct time","recoct size","recoct break","recbip time","recbip size","recbip break"]
        results = DictWriter(f, header)
        results.writeheader()

        # append graph_files
        graph_files = []
        if directory:
            if os.access(filepath, os.W_OK):
                for filename in os.listdir(filepath):
                    if filename.endswith('.s6'):
                        graph_files.append(filename)
        else:
            if '/' not in filepath:
                graph_files.append(filepath)
                filepath = ''
            else:
                graph_files.append(filepath[filepath.rfind('/')+1:len(filepath)])
                filepath = filepath[0:filepath.rfind('/')+1]

        # iterate over graph_files
        for filename in graph_files:
            res = {}

            graphname = filename[0:len(filename)-3]
            print(graphname)
            res["name"] = graphname

            start = time()
            if is_s6(filename):
                graph = read_sparse6("{}{}".format(filepath, filename))
            end = time()
            print("n: {}".format(len(graph)))
            print("time: {}".format(round(end - start, 4)))
            res["n"] = len(graph)

            n = 1
            # approx sol to problem
            if config_args.solve is not None:
                if config_args.solve == 'heuristic':
                    t, minsol, maxsol = run_apx(heuristic_apx, graph, n)
                    print("heuristic apx")
                    print("\tavg time: {}".format(t))
                    print("\tmin size: {}".format(minsol))
                    print("\tmax size: {}".format(maxsol))
                    res["heuristic time"] = t
                    res["heuristic size"] = minsol
                elif config_args.solve == 'dfs':
                    t, minsol, maxsol = run_apx(dfs_apx, graph, n)
                    print("dfs apx")
                    print("\tavg time: {}".format(t))
                    print("\tmin size: {}".format(minsol))
                    print("\tmax size: {}".format(maxsol))
                    res["dfs time"] = t
                    res["dfs size"] = minsol
                elif config_args.solve == 'std':
                    t, minsol, maxsol = run_apx(std_apx, graph, n)
                    print("std apx")
                    print("\tavg time: {}".format(t))
                    print("\tmin size: {}".format(minsol))
                    res["std time"] = t
                    res["std size"] = minsol

            # Edit algorithms
            if config_args.edit == 'remove_octset':
                    start = time()
                    left, right, octset = find_octset(graph)

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
                    res["bip time"] = round(end - start, 4)
                    res["partial"] = len(partial)

                    # lift algorithm
                    if config_args.lift is not None:
                        if config_args.lift == 'greedy':
                            t, minsol, maxsol = run_lift(greedy_lift, graph, n, octset, partial)
                            print("greedy lift")
                            print("\tavg time: {}".format(t))
                            print("\tmin size: {}".format(minsol))
                            print("\tmax size: {}".format(maxsol))
                            res["greedy time"] = t
                            res["greedy size"] = minsol
                        elif config_args.lift == 'naive':
                            t, minsol, maxsol = run_lift(naive_lift, graph, n, octset, partial)
                            print("naive lift")
                            print("\tavg time: {}".format(t))
                            print("\tmin size: {}".format(minsol))
                            print("\tmax size: {}".format(maxsol))
                            res["naive time"] = t
                            res["naive size"] = minsol

            results.writerow(res)
            del graph
            print()


def usage_error_exit(parser, argument, choice=None, list=[]):
    res = 'usage: ' + parser.prog + '\n' + parser.prog + ': error: argument ' + argument
    if choice is None:
        print(res + ': expected one argument')
    else:
        print(res + ': invalid choice: \'' + choice + '\' (choose from ' + str(list)[1:-1] + ')')
    exit(1)


def parse_config():
    edits = ['remove_octset']
    lifts = ['greedy', 'naive']
    solvers = ['dfs', 'heuristic', 'std']

    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s',
        description='Structural Rounding - Experimental Hardness',
        epilog='')
    parser.add_argument('-e', '--edit', choices=edits, help='the editing algorithm')
    parser.add_argument('-l', '--lift', choices=lifts, help='the lifting algorithm')
    parser.add_argument('-s', '--solve', choices=solvers, help='an approximation for the problem')
    parser.add_argument('-g', '--graph', help='the graph file/dir')
    parser.add_argument('-c', '--config', nargs='?', const='config.yaml', help='the optional config (.yaml) file')
    parser.add_argument('-r', '--results', nargs='?', const='results.csv', help='the results (.csv) file')
    parser.add_argument('-v', '--version', action='version', version='Structural Rounding - Experimental Hardness Alpha v1.0')
    config_args = parser.parse_args()

    if config_args.config is not None:
        if is_yaml(config_args.config):
            stream = open(config_args.config, 'r')
            config = yaml.safe_load(stream)

            if config_args.edit is None: config_args.edit = config.get('edit')
            if config_args.lift is None: config_args.lift = config.get('lift')
            if config_args.solve is None: config_args.solve = config.get('solve')
            if config_args.graph is None: config_args.graph = config.get('graph')
            if config_args.results is None: config_args.results = config.get('results')
        else:
            print('error: ' + config_args.config + ' isn\'t a correctly formed (.yaml) file')
            exit(1)

    # check that required arguments exist
    if config_args.graph is None:
        usage_error_exit(parser, '-g/--graph')
    if config_args.results is None:
        usage_error_exit(parser, '-r/--results')

    # check that optional arguments (if from config) are within their choices
    if config_args.solve not in solvers:
        usage_error_exit(parser, '-s/--solve', config_args.solve, solvers)
    if config_args.edit not in edits:
        usage_error_exit(parser, '-e/--edit', config_args.edit, edits)
    if config_args.lift not in lifts:
        usage_error_exit(parser, '-l/--lift', config_args.lift, lifts)

    return config_args


if __name__ == "__main__":
    # cProfile.run("main()")
    main()
