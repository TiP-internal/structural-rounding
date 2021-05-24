
CC=g++
CCFLAGS=-O3 -std=c++11 -fPIC -I./

PYINCLUDE=$(shell python3-config --includes)
PYFLAGS=$(shell python3-config --ldflags) -L. -L./sr_apx/setmap -L./sr_apx/graph -Wl,-rpath,. -Wl,-rpath,./sr_apx/setmap -Wl,-rpath,./sr_apx/graph

# cpp #####################################################################################################

build/util.o: sr_apx/util/util.cpp sr_apx/util/util.hpp
	$(CC) $(CCFLAGS) -c -o build/util.o sr_apx/util/util.cpp

build/matching.o: sr_apx/misc/matching.cpp sr_apx/misc/matching.hpp
	$(CC) $(CCFLAGS) -c -o build/matching.o sr_apx/misc/matching.cpp

build/graph.o: sr_apx/graph/graph.cpp sr_apx/graph/graph.hpp
	$(CC) $(CCFLAGS) -c -o build/graph.o sr_apx/graph/graph.cpp

build/vc_apx.o: sr_apx/vc/apx/vc_apx.cpp sr_apx/vc/apx/vc_apx.hpp
	$(CC) $(CCFLAGS) -c -o build/vc_apx.o sr_apx/vc/apx/vc_apx.cpp

build/vc_exact.o: sr_apx/vc/exact/vc_exact.cpp sr_apx/vc/exact/vc_exact.hpp
	$(CC) $(CCFLAGS) -c -o build/vc_exact.o sr_apx/vc/exact/vc_exact.cpp

build/vc_lift.o: sr_apx/vc/lift/vc_lift.cpp sr_apx/vc/lift/vc_lift.hpp
	$(CC) $(CCFLAGS) -c -o build/vc_lift.o sr_apx/vc/lift/vc_lift.cpp

build/vc_kernel.o: sr_apx/vc/kernel/lp_kernel.cpp sr_apx/vc/kernel/lp_kernel.hpp
	$(CC) $(CCFLAGS) -c -o build/vc_kernel.o sr_apx/vc/kernel/lp_kernel.cpp

build/bipartite.o: sr_apx/bipartite/bipartite.cpp sr_apx/bipartite/bipartite.hpp
	$(CC) $(CCFLAGS) -c -o build/bipartite.o sr_apx/bipartite/bipartite.cpp

build/treewidth.o: sr_apx/treewidth/treewidth.cpp sr_apx/treewidth/treewidth.hpp
	$(CC) $(CCFLAGS) -c -o build/treewidth.o sr_apx/treewidth/treewidth.cpp

build/domset_apx.o: sr_apx/domset/apx/domset_apx.cpp sr_apx/domset/apx/domset_apx.hpp
	$(CC) $(CCFLAGS) -c -o build/domset_apx.o sr_apx/domset/apx/domset_apx.cpp

build/domset_lift.o: sr_apx/domset/lift/domset_lift.cpp sr_apx/domset/lift/domset_lift.hpp
	$(CC) $(CCFLAGS) -c -o build/domset_lift.o sr_apx/domset/lift/domset_lift.cpp

build/domset_exact.o: sr_apx/domset/exact/domset_exact.cpp sr_apx/domset/exact/domset_exact.hpp
	$(CC) $(CCFLAGS) -c -o build/domset_exact.o sr_apx/domset/exact/domset_exact.cpp

lib_sr_apx.so: build/util.o build/matching.o build/graph.o build/vc_apx.o build/vc_exact.o build/vc_lift.o build/vc_kernel.o build/bipartite.o build/treewidth.o build/domset_apx.o build/domset_lift.o build/domset_exact.o sr_apx/setmap/setmap.hpp
	$(CC) -shared -o lib_sr_apx.so build/util.o build/matching.o build/graph.o build/vc_apx.o build/vc_exact.o build/vc_lift.o build/vc_kernel.o build/bipartite.o build/treewidth.o build/domset_apx.o build/domset_lift.o build/domset_exact.o

build/main.o: main.cpp
	$(CC) -O3 -std=c++11 -I./ -c -o build/main.o main.cpp

cpp: build/main.o lib_sr_apx.so
	$(CC) -o main -L. -Wl,-rpath,. build/main.o -l_sr_apx


build/exp.o: exp.cpp
	$(CC) -O3 -std=c++11 -I./ -c -o build/exp.o exp.cpp

exp: build/exp.o lib_sr_apx.so
	$(CC) -o exp -L. -Wl,-rpath,. build/exp.o -l_sr_apx

build/exp1.o: exp1.cpp
	$(CC) -O3 -std=c++11 -I./ -c -o build/exp1.o exp1.cpp

exp1: build/exp1.o lib_sr_apx.so
	$(CC) -o exp1 -L. -Wl,-rpath,. build/exp1.o -l_sr_apx

build/exp2.o: exp2.cpp
	$(CC) -O3 -std=c++11 -I./ -c -o build/exp2.o exp2.cpp

exp2: build/exp2.o lib_sr_apx.so
	$(CC) -o exp2 -L. -Wl,-rpath,. build/exp2.o -l_sr_apx

build/convert.o: convert.cpp
	$(CC) -O3 -std=c++11 -I./ -c -o build/convert.o convert.cpp

convert: build/convert.o lib_sr_apx.so
	$(CC) -o convert -L. -Wl,-rpath,. build/convert.o -l_sr_apx


# python ###########################################################################################################

build/util_module.o: sr_apx/util/util_module.cpp
	$(CC) $(CCFLAGS) -c $(PYINCLUDE) -o build/util_module.o sr_apx/util/util_module.cpp

sr_apx/util/lib_util.so: lib_sr_apx.so build/util_module.o
	$(CC) -shared -o sr_apx/util/lib_util.so build/util_module.o $(PYFLAGS) -l_sr_apx

build/setmap_module.o: sr_apx/setmap/setmap.hpp sr_apx/setmap/setmap_module.cpp sr_apx/setmap/pyset.hpp
	$(CC) $(CCFLAGS) -c $(PYINCLUDE) -o build/setmap_module.o sr_apx/setmap/setmap_module.cpp

sr_apx/setmap/lib_setmap.so: lib_sr_apx.so build/setmap_module.o
	$(CC) -shared -o sr_apx/setmap/lib_setmap.so build/setmap_module.o $(PYFLAGS) -l_sr_apx

build/graph_module.o: sr_apx/graph/graph_module.cpp sr_apx/graph/pygraph.hpp
	$(CC) $(CCFLAGS) -c $(PYINCLUDE) -o build/graph_module.o sr_apx/graph/graph_module.cpp

sr_apx/graph/lib_graph.so: lib_sr_apx.so sr_apx/setmap/lib_setmap.so build/graph_module.o
	$(CC) -shared -o sr_apx/graph/lib_graph.so build/graph_module.o $(PYFLAGS) -l_sr_apx -l_setmap

build/vc_apx_module.o: sr_apx/vc/apx/vc_apx_module.cpp
	$(CC) $(CCFLAGS) -c $(PYINCLUDE) -o build/vc_apx_module.o sr_apx/vc/apx/vc_apx_module.cpp

sr_apx/vc/apx/lib_vc_apx.so: lib_sr_apx.so build/vc_apx_module.o sr_apx/setmap/lib_setmap.so
	$(CC) -shared -o sr_apx/vc/apx/lib_vc_apx.so build/vc_apx_module.o $(PYFLAGS) -l_sr_apx -l_setmap

build/bip_module.o: sr_apx/bipartite/bip_module.cpp
	$(CC) $(CCFLAGS) -c $(PYINCLUDE) -o build/bip_module.o sr_apx/bipartite/bip_module.cpp

sr_apx/bipartite/lib_bipartite.so: lib_sr_apx.so build/bip_module.o sr_apx/setmap/lib_setmap.so
	$(CC) -shared -o sr_apx/bipartite/lib_bipartite.so build/bip_module.o $(PYFLAGS) -l_sr_apx -l_setmap

build/vc_exact_module.o: sr_apx/vc/exact/vc_exact_module.cpp
	$(CC) $(CCFLAGS) -c $(PYINCLUDE) -o build/vc_exact_module.o sr_apx/vc/exact/vc_exact_module.cpp

sr_apx/vc/exact/lib_vc_exact.so: lib_sr_apx.so sr_apx/setmap/lib_setmap.so build/vc_exact_module.o
	$(CC) -shared -o sr_apx/vc/exact/lib_vc_exact.so build/vc_exact_module.o $(PYFLAGS) -l_sr_apx -l_setmap

build/vc_lift_module.o: sr_apx/vc/lift/vc_lift_module.cpp
	$(CC) $(CCFLAGS) -c $(PYINCLUDE) -o build/vc_lift_module.o sr_apx/vc/lift/vc_lift_module.cpp

sr_apx/vc/lift/lib_vc_lift.so: lib_sr_apx.so sr_apx/setmap/lib_setmap.so build/vc_lift_module.o
	$(CC) -shared -o sr_apx/vc/lift/lib_vc_lift.so build/vc_lift_module.o $(PYFLAGS) -l_sr_apx -l_setmap

build/lp_kernel_module.o: sr_apx/vc/kernel/lp_kernel_module.cpp
	$(CC) $(CCFLAGS) -c $(PYINCLUDE) -o build/lp_kernel_module.o sr_apx/vc/kernel/lp_kernel_module.cpp

sr_apx/vc/kernel/lib_lp_kernel.so: lib_sr_apx.so sr_apx/setmap/lib_setmap.so build/lp_kernel_module.o
	$(CC) -shared -o sr_apx/vc/kernel/lib_lp_kernel.so build/lp_kernel_module.o $(PYFLAGS) -l_sr_apx -l_setmap

python: sr_apx/util/lib_util.so sr_apx/setmap/lib_setmap.so sr_apx/graph/lib_graph.so sr_apx/vc/apx/lib_vc_apx.so sr_apx/bipartite/lib_bipartite.so sr_apx/vc/exact/lib_vc_exact.so sr_apx/vc/lift/lib_vc_lift.so sr_apx/vc/kernel/lib_lp_kernel.so

# generator ##########################################################################################

generator: generator/generator.out

generator/generator.out: generator/OCTgenerator.c
	gcc -O3 -std=gnu11 -o generator/generator.out generator/OCTgenerator.c -lm

# remove compiled files ###############################################################################

clean:
	rm -f -r build
	rm -f main
	rm -f exp
	rm -f exp1
	rm -f convert
	rm -f lib_sr_apx.so
	rm -f generator/generator.out
	rm -f sr_apx/util/lib_util.so
	rm -f sr_apx/setmap/lib_setmap.so
	rm -f sr_apx/graph/lib_graph.so
	rm -f sr_apx/vc/apx/lib_vc_apx.so
	rm -f sr_apx/bipartite/lib_bipartite.so
	rm -f sr_apx/vc/exact/lib_vc_exact.so
	rm -f sr_apx/vc/lift/lib_vc_lift.so
	rm -f sr_apx/vc/kernel/lib_lp_kernel.so

# ensures build directory exists when make is invoked
$(shell mkdir -p build)
