
#include <Python.h>

#include "sr_apx/bipartite/bipartite.hpp"
#include "sr_apx/setmap/pyset.hpp"
#include "sr_apx/graph/pygraph.hpp"

static PyObject* bipartite_verifybip(PyObject* self, PyObject* args) {
	PyObject* g;
	PyObject* s;
	if (!PyArg_ParseTuple(args, "OO", &g, &s)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set* os = ((PySet*) s)->s;

	sr_apx::Set octset;
	sr_apx::Set left;
	sr_apx::Set right;

	std::tie(octset, left, right) = sr_apx::bipartite::verify_bipartite(*graph, *os);

	PyObject* l = make_PySet(std::move(left));
	PyObject* r = make_PySet(std::move(right));
	PyObject* o = make_PySet(std::move(octset));
	return Py_BuildValue("OOO", o, l, r);
}

static PyObject* bipartite_prescribed(PyObject* self, PyObject* args) {
	PyObject* bytes;
	PyObject* g;
	if (!PyArg_ParseTuple(args, "OO&", &g, PyUnicode_FSConverter, &bytes)) {
		return NULL;
	}

	char* s;
	Py_ssize_t len;
	PyBytes_AsStringAndSize(bytes, &s, &len);
	Py_DECREF(bytes);

	sr_apx::Graph* graph = ((PyGraph*) g)->g;

	sr_apx::Set oct = sr_apx::bipartite::prescribed_octset(*graph, s);
	PyObject* o = make_PySet(std::move(oct));
	return o;
}

static PyObject* bipartite_vertexdelete(PyObject* self, PyObject* args) {
	PyObject* g;
	if (!PyArg_ParseTuple(args, "O", &g)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;

	sr_apx::Set oct = sr_apx::bipartite::vertex_delete(*graph);
	PyObject* o = make_PySet(std::move(oct));
	return o;
}

static PyMethodDef bipartite_methods[] = {
	{"verify_bipartite", bipartite_verifybip, METH_VARARGS, "computes an oct decomposition from a given octset"},
	{"prescribed_octset", bipartite_prescribed, METH_VARARGS, "reads a predetermined octset from a file"},
	{"vertex_delete", bipartite_vertexdelete, METH_VARARGS, "computes an octset for a given graph"},
	{NULL},
};

static struct PyModuleDef bipartite_module = {
	PyModuleDef_HEAD_INIT,
	"octset",
	"Python interface for oct finding operations",
	-1,
	bipartite_methods
};

PyMODINIT_FUNC PyInit_lib_bipartite() {
	return PyModule_Create(&bipartite_module);
}
