
#include <Python.h>

#include "sr_apx/vc/exact/vc_exact.hpp"
#include "sr_apx/graph/pygraph.hpp"
#include "sr_apx/setmap/pyset.hpp"

static PyObject* vc_exact_bipexact(PyObject* self, PyObject* args) {
	PyObject* g;
	if (!PyArg_ParseTuple(args, "O", &g)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set cover = sr_apx::vc::exact::bip_exact(*graph);
	return make_PySet(std::move(cover));
}

static PyMethodDef vc_exact_methods[] = {
	{"bip_exact", vc_exact_bipexact, METH_VARARGS, "computes a minimum vertex cover in a bipartite graph"},
	{NULL},
};

static struct PyModuleDef vc_exact_module = {
	PyModuleDef_HEAD_INIT,
	"vc_exact",
	"Python interface for vertex cover algorithms",
	-1,
	vc_exact_methods
};

PyMODINIT_FUNC PyInit_lib_vc_exact() {
	return PyModule_Create(&vc_exact_module);
}
