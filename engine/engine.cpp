#include <pybind11/pybind11.h>

#include "headers/Settings.h"

namespace py = pybind11;

PYBIND11_MODULE(engine, m) {
    py::class_<BondParameters>(m, "BondParameters")
            .def(py::init<>())
            .def_readwrite("lambda_", &BondParameters::lambda_)
            .def_readwrite("sigma", &BondParameters::sigma)
            ;
}