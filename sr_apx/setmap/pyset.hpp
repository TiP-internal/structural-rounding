
#ifndef PYSET_H
#define PYSET_H

#include <Python.h>
#include "sr_apx/setmap/setmap.hpp"

typedef struct {
	PyObject_HEAD
	Set* s;
	bool borrowed;
} PySet;

PyObject* make_PySet(Set*, bool);

#endif
