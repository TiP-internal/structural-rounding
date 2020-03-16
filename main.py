# header = ["name","n","m","dfs time","dfs size","heuristic time","heuristic size","std time","std size","stdrev time","stdrev size","oct size","partial","bip time","naive time","naive size","apx time","apx size","greedy time","greedy size","octfirst time","octfirst size","octfirst break","bipfirst time","bipfirst size","bipfirst break","rec time","rec size","rec break","recoct time","recoct size","recoct break","recbip time","recbip size","recbip break"]

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


def print_result(algo, time, min_sol, max_sol):
    print(algo)
    print("\tavg time: {}\n\tmin size: {}\n\tmax size: {}"
        .format(time, min_sol, max_sol))


def main():
    # parse config
    pipes = parse_config()
    print(pipes.values())
    exit()

    # iterate over pipes
    first = True
    for pipe in pipes.values():
        if first:
            first = False
        else:
            print('\n')

        # format filepath
        filepath = pipe['graph']
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

        # write to results file
        with open(pipe['results'], "w") as f:
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

            header = ["name","n","m"]

            # iterate over graph_files
            first = True
            for filename in graph_files:
                if first:
                    first = False
                else:
                    print()

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

                # approximate the problem
                n = 1
                solve_algo = pipe['solve']
                solve_time, solve_size = solve_algo + ' time', solve_algo + ' size'
                header.extend([solve_time, solve_size])

                if pipe['solve'] == 'heuristic':
                    t, minsol, maxsol = run_apx(heuristic_apx, graph, n)
                elif pipe['solve'] == 'dfs':
                    t, minsol, maxsol = run_apx(dfs_apx, graph, n)
                elif pipe['solve'] == 'std':
                    t, minsol, maxsol = run_apx(std_apx, graph, n)

                print_result(solve_algo + ' apx', t, minsol, maxsol)
                res[solve_time] = t
                res[solve_size] = minsol

                # edit
                if pipe['edit'] == 'remove_octset':
                    start = time()
                    left, right, octset = find_octset(graph)

                    bippart = Set()
                    for v in left:
                        bippart.add(v)
                    for v in right:
                        bippart.add(v)

                    # SR solve
                    partial = bip_exact(graph.subgraph(bippart))
                    end = time()

                    print("bip solve")
                    print("\tavg time: {}".format(round(end - start, 4)))
                    print(len(partial))

                    header.extend(['bip time', 'oct size', 'partial'])
                    res["oct size"] = len(octset)
                    res["bip time"] = round(end - start, 4)
                    res["partial"] = len(partial)

                    # lift
                    if pipe['lift'] is not None:
                        lift_algo = pipe['lift'] # We may assume up to this point `lift` has been defined
                        lift_time, lift_size = lift_algo + ' time', lift_algo + ' size'
                        header.extend([lift_time, lift_size])

                        if pipe['lift'] == 'greedy':
                            t, minsol, maxsol = run_lift(greedy_lift, graph, n, octset, partial)
                        elif pipe['lift'] == 'naive':
                            t, minsol, maxsol = run_lift(naive_lift, graph, n, octset, partial)

                        print_result(lift_algo + ' lift', t, minsol, maxsol)
                        res[lift_time] = t
                        res[lift_size] = minsol

                results = DictWriter(f, header)
                results.writeheader()

                results.writerow(res)
                del graph


def usage_error_exit(parser, argument, choice=None, list=[]):
    res = 'usage: {}\n{}: error: argument {}'.format(parser.prog, parser.prog, argument)
    if choice is None:
        print(res + ': expected one argument')
    else:
        print(res + ': invalid choice: \'{}\' (choose from {})'.format(choice, str(list)[1:-1]))
    exit(1)


def parse_config():
    edits   = ['remove_octset']
    solvers = ['dfs', 'heuristic', 'std']
    lifts   = ['greedy', 'naive']

    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s',
        description='Structural Rounding - Experimental Harness',
        epilog='')
    parser.add_argument('-g', '--graph', help='the graph file/dir')
    parser.add_argument('-r', '--results', nargs='?', const='results.csv', help='the results (.csv) file')
    parser.add_argument('-c', '--config', nargs='?', const='config.yaml', help='the optional config (.yaml) file')
    parser.add_argument('-e', '--edit', choices=edits, help='the editing algorithm')
    parser.add_argument('-s', '--solve', choices=solvers, help='an approximation for the problem')
    parser.add_argument('-l', '--lift', choices=lifts, help='the lifting algorithm')
    parser.add_argument('-v', '--version', action='version', version='Structural Rounding - Experimental Harness, Alpha v1.0')
    config_args = parser.parse_args()

    if config_args.config is not None:
        if is_yaml(config_args.config):
            pipes = {}

            # single-pipe config file
            try:
                pipe = yaml.safe_load(open(config_args.config, 'r'))
                if config_args.graph is not None: pipe['graph'] = config_args.graph
                if config_args.results is not None: pipe['results'] = config_args.results
                if config_args.edit is not None: pipe['edit'] = config_args.edit
                if config_args.solve is not None: pipe['solve'] = config_args.solve
                if config_args.lift is not None: pipe['lift'] = config_args.lift
                pipes[0] = pipe

            # multi-pipe config file
            except:
                i = 0;
                for pipe in yaml.safe_load_all(open(config_args.config, 'r')):
                    if config_args.graph is not None: pipe['graph'] = config_args.graph
                    if config_args.results is not None: pipe['results'] = config_args.results
                    if config_args.edit is not None: pipe['edit'] = config_args.edit
                    if config_args.solve is not None: pipe['solve'] = config_args.solve
                    if config_args.lift is not None: pipe['lift'] = config_args.lift
                    pipes[i] = pipe
                    i += 1

        else:
            print('error: ' + config_args.config + ' isn\'t a correctly formed (.yaml) file')
            exit(1)
    else:
        pipes = {}
        pipes[0] = {'graph': config_args.graph, 'edit': config_args.edit, 'solve': config_args.solve, 'lift': config_args.lift, 'results': config_args.results}

    for pipe in pipes.values():
        # check that required arguments exist
        if 'graph' not in pipe.keys() or pipe['graph'] is None:
            usage_error_exit(parser, '-g/--graph')
        if 'results' not in pipe.keys() or pipe['results'] is None:
            usage_error_exit(parser, '-r/--results')

        # check that optional arguments (if from config) are within their choices
        if ('edit' in pipe.keys() and pipe['edit'] is not None) and (pipe['edit'] not in edits):
            usage_error_exit(parser, '-e/--edit', pipe['edit'], edits)
        if ('solve' in pipe.keys() and pipe['solve'] is not None) and (pipe['solve'] not in solvers):
            usage_error_exit(parser, '-s/--solve', pipe['solve'], solvers)
        if ('lift' in pipe.keys() and pipe['lift'] is not None) and (pipe['lift'] not in lifts):
            usage_error_exit(parser, '-l/--lift', pipe['lift'], lifts)

    return pipes


if __name__ == "__main__":
    # cProfile.run("main()")
    main()
