
#ifndef PYSET_H
#define PYSET_H

#include <Python.h>

#include "sr_apx/setmap/setmap.hpp"

typedef struct {
	PyObject_HEAD
	sr_apx::Set* s;
	bool borrowed;
} PySet;

PyObject* make_PySet(sr_apx::Set&&);
PyObject* make_PySet(sr_apx::Set*);

#endif
