
#ifndef PYGRAPH_H
#define PYGRAPH_H

#include <Python.h>

#include "sr_apx/graph/graph.hpp"

typedef struct {
	PyObject_HEAD
	sr_apx::Graph* g;
} PyGraph;

PyObject* make_PyGraph(sr_apx::Graph&&);

#endif
