# header = ["name","n","m","dfs time","dfs size","heuristic time","heuristic size","std time","std size","stdrev time","stdrev size","oct size","partial","bip time","naive time","naive size","apx time","apx size","greedy time","greedy size","octfirst time","octfirst size","octfirst break","bipfirst time","bipfirst size","bipfirst break","rec time","rec size","rec break","recoct time","recoct size","recoct break","recbip time","recbip size","recbip break"]

import cProfile

import sys, os, argparse, yaml
from time import time
from csv import DictWriter
from multiprocessing import Process

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


def process_files(pipe):
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

    # append graph files
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

    return filepath, graph_files


def graph_dependencies(pipes):
    approximate_solves = {'heuristic', 'dfs', 'std'}
    graph_pipes = {}

    for pipe in pipes.values():
        filepath, graph_files = process_files(pipe)
        for file in graph_files:
            fullpath = filepath+file

            # approximate solves
            if pipe['solve'] in approximate_solves:
                if fullpath+pipe['solve'] in graph_pipes.keys():
                    graph_pipes[fullpath+pipe['solve']] = 1
                else:
                    graph_pipes[fullpath+pipe['solve']] = None

            # exact solves
            else:
                if fullpath+pipe['edit'] in graph_pipes.keys():
                    graph_pipes[fullpath+pipe['edit']] = 1
                    if fullpath+pipe['edit']+pipe['solve'] in graph_pipes.keys():
                        graph_pipes[fullpath+pipe['edit']+pipe['solve']] = 1
                        if fullpath+pipe['edit']+pipe['solve']+pipe['lift'] in graph_pipes.keys():
                            graph_pipes[fullpath+pipe['edit']+pipe['solve']+pipe['lift']] = 1
                        else:
                            graph_pipes[fullpath+pipe['edit']+pipe['solve']+pipe['lift']] = None
                    else:
                        graph_pipes[fullpath+pipe['edit']+pipe['solve']] = None
                        graph_pipes[fullpath+pipe['edit']+pipe['solve']+pipe['lift']] = None
                else:
                    graph_pipes[fullpath+pipe['edit']] = None
                    graph_pipes[fullpath+pipe['edit']+pipe['solve']] = None
                    graph_pipes[fullpath+pipe['edit']+pipe['solve']+pipe['lift']] = None

    del_list = set()
    for pipe in graph_pipes.keys():
        if graph_pipes[pipe] == None:
            del_list.add(pipe)
    for unique in del_list:
        del graph_pipes[unique]

    return graph_pipes


def main():
    # parse config and find redundant work
    pipes = parse_config()
    duplicated_steps = graph_dependencies(pipes)

    # iterate over pipes
    first = True
    for pipe in pipes.values():
        if first: first = False
        else: print('\n')

        # get filepath and graph_files
        filepath, graph_files = process_files(pipe)

        # write to results file
        with open(pipe['results'], "w") as f:
            header = ["name","n","m"]

            # iterate over graph_files
            first = True
            for filename in graph_files:
                if first: first = False
                else: print()

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

                # run solution pipeline on another process
                pipeline = Process(target=run_pipeline, args=(pipe, filepath, filename, duplicated_steps, graph, header, n, res))
                pipeline.start()

                # kill solution pipeline process if timed out
                pipeline.join(timout=pipe['timeout'])

                results = DictWriter(f, header)
                results.writeheader()

                # if timed out
                if pipeline.exitcode is None:
                    # write timout error withe line no information


                    # fill empty places in results with '-'
                    for key in header:
                        if key not in res.keys():
                            res[key] = '-'

                results.writerow(res)
                del graph

def run_pipeline(pipe, filepath, filename, duplicated_steps, graph, header, n, res):
    # edit
    if pipe['edit'] is not None:
        edit_key = filepath + filename + pipe['edit']
        if edit_key in duplicated_steps:
            if duplicated_steps[edit_key] == 1:
                if pipe['edit'] == 'remove_octset':
                    start = time()
                    left, right, octset = find_octset(graph)
                    bippart = Set()
                    for v in left:
                        bippart.add(v)
                    for v in right:
                        bippart.add(v)
                    duplicated_steps[edit_key] = bippart
            else:
                bippart = duplicated_steps[edit_key]
        else:
            if pipe['edit'] == 'remove_octset':
                start = time()
                left, right, octset = find_octset(graph)
                bippart = Set()
                for v in left:
                    bippart.add(v)
                for v in right:
                    bippart.add(v)

    # solve
    if pipe['solve'] is not None:
        solve_algo = pipe['solve']
        solve_time, solve_size = solve_algo + ' time', solve_algo + ' size'
        header.extend([solve_time, solve_size])

        # approximate solve
        approximate_solves = {'heuristic', 'dfs', 'std'}
        if solve_algo in approximate_solves:
            approx_solve_key = filepath + filename + solve_algo
            if approx_solve_key in duplicated_steps:
                if duplicated_steps[approx_solve_key] == 1:
                    if solve_algo == 'heuristic':
                        t, minsol, maxsol = run_apx(heuristic_apx, graph, n)
                    elif solve_algo == 'dfs':
                        t, minsol, maxsol = run_apx(dfs_apx, graph, n)
                    elif solve_algo == 'std':
                        t, minsol, maxsol = run_apx(std_apx, graph, n)
                    duplicated_steps[approx_solve_key] = t, minsol, maxsol
                else:
                    t, minsol, maxsol = duplicated_steps[approx_solve_key]
            else:
                if pipe['solve'] == 'heuristic':
                    t, minsol, maxsol = run_apx(heuristic_apx, graph, n)
                elif pipe['solve'] == 'dfs':
                    t, minsol, maxsol = run_apx(dfs_apx, graph, n)
                elif pipe['solve'] == 'std':
                    t, minsol, maxsol = run_apx(std_apx, graph, n)

            print_result(solve_algo + ' apx', t, minsol, maxsol)
            res[solve_time] = t
            res[solve_size] = minsol

        # exact solve
        else:
            solve_key = edit_key + solve_algo
            if solve_algo == 'bip_exact':
                if solve_key in duplicated_steps:
                    if duplicated_steps[solve_key] == 1:
                        partial = bip_exact(graph.subgraph(bippart))
                        end = time()
                        t = end - start
                        duplicated_steps[solve_key] = partial, t
                    else:
                        partial, t = duplicated_steps[solve_key]
                else:
                    partial = bip_exact(graph.subgraph(bippart))
                    end = time()
                    t = end - start

                print("bip solve")
                print("\tavg time: {}".format(round(t, 4)))
                print(len(partial))

                header.extend(['bip time', 'oct size', 'partial'])
                res["oct size"] = len(octset)
                res["bip time"] = round(t, 4)
                res["partial"] = len(partial)

    # lift
    if pipe['lift'] is not None and pipe['solve'] not in approximate_solves:
        lift_algo = pipe['lift']  # We may assume up to this point `lift` has been defined
        lift_key = solve_key + lift_algo
        lift_time, lift_size = lift_algo + ' time', lift_algo + ' size'
        header.extend([lift_time, lift_size])

        if lift_key in duplicated_steps:
            if duplicated_steps[lift_key] == 1:
                if pipe['lift'] == 'greedy':
                    t, minsol, maxsol = run_lift(greedy_lift, graph, n, octset, partial)
                elif pipe['lift'] == 'naive':
                    t, minsol, maxsol = run_lift(naive_lift, graph, n, octset, partial)
                duplicated_steps[lift_key] = t, minsol, maxsol
            else:
                t, minsol, maxsol = duplicated_steps[lift_key]
        else:
            if pipe['lift'] == 'greedy':
                t, minsol, maxsol = run_lift(greedy_lift, graph, n, octset, partial)
            elif pipe['lift'] == 'naive':
                t, minsol, maxsol = run_lift(naive_lift, graph, n, octset, partial)

        print_result(lift_algo + ' lift', t, minsol, maxsol)
        res[lift_time] = t
        res[lift_size] = minsol

def usage_error_exit(parser, argument, choice=None, list=[]):
    res = 'usage: {}\n{}: error: argument {}'.format(parser.prog, parser.prog, argument)
    if choice is None:
        print(res + ': expected one argument')
    else:
        print(res + ': invalid choice: \'{}\' (choose from {})'.format(choice, str(list)[1:-1]))
    exit(1)


def parse_config():
    edits   = ['remove_octset']
    solvers = ['dfs', 'heuristic', 'std', 'bip_exact']
    lifts   = ['greedy', 'naive']

    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s',
        description='Structural Rounding - Experimental Harness',
        epilog='')
    parser.add_argument('-g', '--graph', help='the graph file/dir')
    parser.add_argument('-c', '--config', nargs='?', const='config.yaml', help='the optional config (.yaml) file')
    parser.add_argument('-r', '--results', nargs='?', const='results.csv', help='the results (.csv) file')
    parser.add_argument('-e', '--edit', choices=edits, help='the editing algorithm')
    parser.add_argument('-s', '--solve', choices=solvers, help='an approximation/or solution for the problem')
    parser.add_argument('-l', '--lift', choices=lifts, help='the lifting algorithm')
    parser.add_argument('-v', '--version', action='version', version='Structural Rounding - Experimental Harness, Alpha v1.0')
    parser.add_argument('-t', '--timeout', default=3600, help='sets the timeout, in seconds, for each pipeline')
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
                if config_args.timeout is not None: pipe['timeout'] = config_args.timeout
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
                    if config_args.timeout is not None: pipe['timeout'] = config_args.timeout
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
