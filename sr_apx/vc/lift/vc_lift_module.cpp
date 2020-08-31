
#include <Python.h>

#include "sr_apx/vc/lift/vc_lift.hpp"
#include "sr_apx/graph/pygraph.hpp"
#include "sr_apx/setmap/pyset.hpp"

static PyObject* vc_lift_naivelift(PyObject* self, PyObject* args) {
	PyObject* g;
	PyObject* o;
	PyObject* p;

	if (!PyArg_ParseTuple(args, "OOO", &g, &o, &p)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set* octset = ((PySet*) o)->s;
	sr_apx::Set* partial = ((PySet*) p)->s;

	return make_PySet(sr_apx::vc::lift::naive_lift(graph, octset, partial), false);
}

static PyObject* vc_lift_greedylift(PyObject* self, PyObject* args) {
	PyObject* g;
	PyObject* o;
	PyObject* p;

	if (!PyArg_ParseTuple(args, "OOO", &g, &o, &p)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set* octset = ((PySet*) o)->s;
	sr_apx::Set* partial = ((PySet*) p)->s;

	return make_PySet(sr_apx::vc::lift::greedy_lift(graph, octset, partial), false);
}

static PyObject* vc_lift_apxlift(PyObject* self, PyObject* args) {
	PyObject* g;
	PyObject* o;
	PyObject* p;

	if (!PyArg_ParseTuple(args, "OOO", &g, &o, &p)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set* octset = ((PySet*) o)->s;
	sr_apx::Set* partial = ((PySet*) p)->s;

	return make_PySet(sr_apx::vc::lift::apx_lift(graph, octset, partial), false);
}

static PyObject* vc_lift_octlift(PyObject* self, PyObject* args) {
	PyObject* g;
	PyObject* o;
	PyObject* p;

	if (!PyArg_ParseTuple(args, "OOO", &g, &o, &p)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set* octset = ((PySet*) o)->s;
	sr_apx::Set* partial = ((PySet*) p)->s;

	return make_PySet(sr_apx::vc::lift::oct_lift(graph, octset, partial), false);
}

static PyObject* vc_lift_biplift(PyObject* self, PyObject* args) {
	PyObject* g;
	PyObject* o;
	PyObject* p;

	if (!PyArg_ParseTuple(args, "OOO", &g, &o, &p)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set* octset = ((PySet*) o)->s;
	sr_apx::Set* partial = ((PySet*) p)->s;

	return make_PySet(sr_apx::vc::lift::bip_lift(graph, octset, partial), false);
}

static PyObject* vc_lift_recursivelift(PyObject* self, PyObject* args) {
	PyObject* g;
	PyObject* o;
	PyObject* p;

	if (!PyArg_ParseTuple(args, "OOO", &g, &o, &p)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set* octset = ((PySet*) o)->s;
	sr_apx::Set* partial = ((PySet*) p)->s;

	return make_PySet(sr_apx::vc::lift::recursive_lift(graph, octset, partial), false);
}

static PyObject* vc_lift_recoctlift(PyObject* self, PyObject* args) {
	PyObject* g;
	PyObject* o;
	PyObject* p;

	if (!PyArg_ParseTuple(args, "OOO", &g, &o, &p)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set* octset = ((PySet*) o)->s;
	sr_apx::Set* partial = ((PySet*) p)->s;

	return make_PySet(sr_apx::vc::lift::recursive_oct_lift(graph, octset, partial), false);
}

static PyObject* vc_lift_recbiplift(PyObject* self, PyObject* args) {
	PyObject* g;
	PyObject* o;
	PyObject* p;

	if (!PyArg_ParseTuple(args, "OOO", &g, &o, &p)) {
		return NULL;
	}

	sr_apx::Graph* graph = ((PyGraph*) g)->g;
	sr_apx::Set* octset = ((PySet*) o)->s;
	sr_apx::Set* partial = ((PySet*) p)->s;

	return make_PySet(sr_apx::vc::lift::recursive_bip_lift(graph, octset, partial), false);
}

static PyMethodDef vc_lift_methods[] = {
	{"naive_lift", vc_lift_naivelift, METH_VARARGS, "computes a complete vertex cover from a partial solution and an octset"},
	{"greedy_lift", vc_lift_greedylift, METH_VARARGS, "computes a complete vertex cover from a partial solution and an octset"},
	{"apx_lift", vc_lift_apxlift, METH_VARARGS, "computes a complete vertex cover from a partial solution and an octset"},
	{"oct_lift", vc_lift_octlift, METH_VARARGS, "computes a complete vertex cover from a partial solution and an octset"},
	{"bip_lift", vc_lift_biplift, METH_VARARGS, "computes a complete vertex cover from a partial solution and an octset"},
	{"recursive_lift", vc_lift_recursivelift, METH_VARARGS, "computes a complete vertex cover from a partial solution and an octset"},
	{"recursive_oct_lift", vc_lift_recoctlift, METH_VARARGS, "computes a complete vertex cover from a partial solution and an octset"},
	{"recursive_bip_lift", vc_lift_recbiplift, METH_VARARGS, "computes a complete vertex cover from a partial solution and an octset"},
	{NULL},
};

static struct PyModuleDef vc_lift_module = {
	PyModuleDef_HEAD_INIT,
	"vc_lift",
	"Python interface for vertex cover algorithms",
	-1,
	vc_lift_methods
};

PyMODINIT_FUNC PyInit_lib_vc_lift() {
	return PyModule_Create(&vc_lift_module);
}
