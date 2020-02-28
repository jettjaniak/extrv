#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// for string formatting
#include <iomanip>
#include <sstream>

#include "types.h"
#include "Settings.h"
#include "SimulationState.h"

namespace py = pybind11;
using pybind11::literals::operator""_a;

// TODO: __repr__ for each object

PYBIND11_MODULE(extrv_engine, m) {

    /////////////
    //  types  //
    /////////////

    py::class_<xy_t> xy_type(m, "xy_t");

    xy_type.def("__repr__",
        [](const xy_t &xy) {
            std::stringstream stream;
            stream << std::fixed << std::setprecision(4);
            stream << "<xy_t (" << xy.x << ", " << xy.y << ")>";
            return stream.str();
        }
    );

    xy_type.def_readwrite("x", &xy_t::x)
           .def_readwrite("y", &xy_t::y);


    //////////////////////
    //  BondParameters  //
    //////////////////////

    py::class_<BondParameters> bond_p(m, "BondParameters");

    py::enum_<BondParameters::BondType >(bond_p, "BondType")
            .value("psel", BondParameters::BondType::psel)
            .value("esel", BondParameters::BondType::esel)
            .value("integrin", BondParameters::BondType::integrin);

    bond_p.def(
        py::init<BondParameters::BondType, double, double, double, double, double, double, double, double>(),
        "bond_type"_a, "lambda_"_a, "sigma"_a, "k_f_0"_a, "rec_dens"_a, "x1s"_a, "k01s"_a, "x1c"_a=0.0, "k01c"_a=0.0
    );

    bond_p.def_readwrite("bond_type", &BondParameters::bond_type)
          .def_readwrite("lambda_", &BondParameters::lambda_)
          .def_readwrite("sigma", &BondParameters::sigma)
          .def_readwrite("k_f_0", &BondParameters::k_f_0)
          .def_readwrite("rec_dens", &BondParameters::rec_dens)
          .def_readwrite("x1s", &BondParameters::x1s)
          .def_readwrite("k01s", &BondParameters::k01s)
          .def_readwrite("x1c", &BondParameters::x1c)
          .def_readwrite("k01c", &BondParameters::k01c);


    //////////////////
    //  LigandType  //
    //////////////////

    py::class_<LigandType> lig_type(m, "LigandType");

    // passing a pointer is OK
    lig_type.def(py::init())
            .def("add_bond_p", &LigandType::add_bond_p, "bond_p"_a);

    lig_type.def_readonly("bonds_p", &LigandType::bonds_p);


    ///////////////////////
    //  ModelParameters  //
    ///////////////////////

    py::class_<ModelParameters> p(m, "ModelParameters");

    p.def(
        py::init<double, double, double, double, double, double>(),
        "r_c"_a, "mu"_a, "temp"_a, "dens_diff"_a, "f_rep_0"_a, "tau"_a
    );

    p.def_readwrite("r_c", &ModelParameters::r_c)
     .def_readwrite("mu", &ModelParameters::mu)
     .def_readwrite("temp", &ModelParameters::temp)
     .def_readwrite("dens_diff", &ModelParameters::dens_diff)
     .def_readwrite("f_rep_0", &ModelParameters::f_rep_0)
     .def_readwrite("tau", &ModelParameters::tau);


    ////////////////
    //  Settings  //
    ////////////////

    py::class_<Settings> settings(m, "Settings");

    settings.def(py::init<ModelParameters*>(), "p"_a)
            .def("add_lig_type", &Settings::add_lig_type, "lig_type"_a, "n_of_lig"_a);

    settings.def_readwrite("p", &Settings::p)
            .def_readwrite("lig_types_and_nrs", &Settings::lig_types_and_nrs);


    ///////////////
    //  History  //
    ///////////////

    py::class_<History> hist(m, "History");

    py::class_<History::BondTrajectory> bond_traj(hist, "BondTrajectory");

    bond_traj.def_readonly("start_i", &History::BondTrajectory::start_i)
            .def_readonly("positions", &History::BondTrajectory::positions);

    hist.def_readonly("h", &History::h)
            .def_readonly("rot", &History::rot)
            .def_readonly("bond_trajectories", &History::bond_trajectories);


    ///////////////////////
    //  SimulationState  //
    ///////////////////////

    py::class_<SimulationState> s(m, "SimulationState");

    s.def(
        py::init<double, Settings*, unsigned int>(),
        "h_0"_a, "settings"_a, "seed"_a
    );

    s.def("simulate_one_step", &SimulationState::simulate_one_step, "dt"_a, "shear"_a)
     .def("simulate", &SimulationState::simulate, "n_steps"_a, "dt"_a, "shear"_a)
     .def("simulate_with_history", &SimulationState::simulate_with_history,
          "n_steps"_a, "dt"_a, "shear"_a, "save_every"_a=1000);

    s.def_readwrite("h", &SimulationState::h)
     .def_readwrite("rot", &SimulationState::rot)
     .def_readwrite("bd_lig_ind", &SimulationState::bd_lig_ind)
     .def_readwrite("settings", &SimulationState::settings);



}


