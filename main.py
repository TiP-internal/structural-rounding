
import cProfile

import os, argparse, yaml
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


def main():
    # parse config
    config_args = parse_config()

    # format filepath
    filepath = config_args.graph
    directory = True
    if filepath[len(filepath)-1] != '/':
        if is_s6(filepath):
            directory = False
        elif os.path.isdir(os.path.join(os.getcwd(), filepath)):
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

            # problem
            if config_args.problem == 'vertex_cover':
                n = 1
                # approx sol to problem
                if config_args.approx is not None:
                    if config_args.approx == 'heuristic':
                        t, minsol, maxsol = run_apx(heuristic_apx, graph, n)
                        print("heuristic apx")
                        print("\tavg time: {}".format(t))
                        print("\tmin size: {}".format(minsol))
                        print("\tmax size: {}".format(maxsol))
                        res["heuristic time"] = t
                        res["heuristic size"] = minsol
                    elif config_args.approx == 'dfs':
                        t, minsol, maxsol = run_apx(dfs_apx, graph, n)
                        print("dfs apx")
                        print("\tavg time: {}".format(t))
                        print("\tmin size: {}".format(minsol))
                        print("\tmax size: {}".format(maxsol))
                        res["dfs time"] = t
                        res["dfs size"] = minsol
                    elif config_args.approx == 'std':
                        t, minsol, maxsol = run_apx(std_apx, graph, n)
                        print("std apx")
                        print("\tavg time: {}".format(t))
                        print("\tmin size: {}".format(minsol))
                        res["std time"] = t
                        res["std size"] = minsol

                # class
                if config_args.gclass == 'bipartite':
                    # edit algorithm
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


def usage_error(parser, argument, info, choices=None):
    if choices is None:
        return parser.prog + ': error: argument ' + argument + ': ' + info
    else:
        return parser.prog + ': error: argument ' + argument + ': ' + info + ' (choose from ' + choices + ')'


def parse_config():
    problems = ['vertex_cover']
    classes = ['bipartite']
    edits = ['remove_octset']
    lifts = ['greedy', 'naive']
    approx = ['dfs', 'heuristic', 'std']

    parser = argparse.ArgumentParser(prog='main.py', usage='%(prog)s',
        description='Structural Rounding - Experimental Harness',
        epilog='')
    parser.add_argument('-p', '--problem', choices=problems, help='the problem')
    parser.add_argument('-c', '--class', dest='gclass', choices=classes, help='the graph class to edit to')
    parser.add_argument('-e', '--edit', choices=edits, help='the editing algorithm')
    parser.add_argument('-l', '--lift', choices=lifts, help='the lifting algorithm')
    parser.add_argument('-a', '--approx', choices=approx, help='an approximation for the problem')
    parser.add_argument('-g', '--graph', help='the graph file/dir')
    parser.add_argument('-s', '--spec', nargs='?', const='config.yaml', help='the optional config (.yaml) file')
    parser.add_argument('-r', '--results', nargs='?', const='results.csv', help='the results (.csv) file')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s Alpha v1.0')
    config_args = parser.parse_args()

    if config_args.spec is not None:
        if is_yaml(config_args.spec):
            stream = open(config_args.spec, 'r')
            config = yaml.safe_load(stream)
            if config_args.problem is None: config_args.problem = config.get('problem')
            if config_args.gclass is None: config_args.gclass = config.get('class')
            if config_args.edit is None: config_args.edit = config.get('edit')
            if config_args.lift is None: config_args.lift = config.get('lift')
            if config_args.approx is None: config_args.approx = config.get('approx')
            if config_args.graph is None: config_args.graph = config.get('graph')
            if config_args.results is None: config_args.results = config.get('results')
        else:
            print('error: ' + config_args.spec + ' isn\'t a correctly formed (.yaml) file')
            exit(1)

    # check that required arguments exist
    if config_args.problem is None:
        print(usage_error(parser, '-p/--problem', 'expected one argument'))
        exit(1)
    if config_args.problem is not None and config_args.problem != 'vertex_cover':
        # print(usage_error(parser, '-p/--problem', 'invalid choice', problems))
        # print('usage: ' + parser.prog + '\n' + parser.prog + ': error: argument -p/--problem: invalid choice: \'' + config_args.problem + '\' (choose from \'vertex_cover\')')
        exit(1)


    if config_args.graph is None:
        print('usage: ' + parser.prog + '\n' + parser.prog + ': error: argument -g/--graph: expected one argument')
        exit(1)
    if config_args.results is None:
        print('usage: ' + parser.prog + '\n' + parser.prog + ': error: argument -r/--results: expected one argument')
        exit(1)

    # check that optional arguments are within their choices
    if config_args.problem is not None and config_args.problem != 'vertex_cover':
        # print(usage_error(parser, '-p/--problem', 'invalid choice', problems))
        # print('usage: ' + parser.prog + '\n' + parser.prog + ': error: argument -p/--problem: invalid choice: \'' + config_args.problem + '\' (choose from \'vertex_cover\')')
        exit(1)
    if config_args.gclass is not None and config_args.gclass != 'bipartite':
        print('usage: ' + parser.prog + '\n' + parser.prog + ': error: argument -c/--class: invalid choice: \'' + config_args.gclass + '\' (choose from \'bipartite\')')
        exit(1)
    if config_args.edit is not None and config_args.edit != 'remove_octset':
        print('usage: ' + parser.prog + '\n' + parser.prog + ': error: argument -e/--edit: invalid choice: \'' + config_args.edit + '\' (choose from \'remove_octset\')')
        exit(1)
    if config_args.lift is not None and config_args.lift != 'greedy' and config_args.lift != 'naive':
        print('usage: ' + parser.prog + '\n' + parser.prog + ': error: argument -l/--lift: invalid choice: \'' + config_args.lift + '\' (choose from \'greedy\', \'naive\')')
        exit(1)
    if config_args.approx is not None and config_args.approx != 'dfs' and config_args.approx != 'heuristic' and config_args.approx != 'std':
        print('usage: ' + parser.prog + '\n' + parser.prog + ': error: argument -a/--approx: invalid choice: \'' + config_args.approx + '\' (choose from \'dfs\', \'heuristic\', \'std\')')
        exit(1)

    return config_args


def is_yaml(file):
    return file.lower().endswith(('.yaml'))


def is_s6(file):
    return file.lower().endswith(('.s6'))


if __name__ == "__main__":
    # cProfile.run("main()")
    main()
