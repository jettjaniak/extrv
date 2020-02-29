#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// for string formatting
#include <iomanip>
#include <sstream>

#include "types.h"
#include "Settings.h"
#include "SimulationState.h"
#include "History.h"
#include "AbstractBondType.h"

namespace py = pybind11;
using pybind11::literals::operator""_a;

// TODO: __repr__ for each object (if it makes sense)

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


    ////////////////////////
    //  AbstractBondType  //
    ////////////////////////

    py::class_<AbstractBondType> abstr_bond_t(m, "AbstractBondType");

    abstr_bond_t.def_readwrite("lambda_", &AbstractBondType::lambda_)
                .def_readwrite("sigma", &AbstractBondType::sigma)
                .def_readwrite("k_f_0", &AbstractBondType::k_f_0)
                .def_readwrite("rec_dens", &AbstractBondType::rec_dens)
                .def_readwrite("x1s", &AbstractBondType::x1s)
                .def_readwrite("k01s", &AbstractBondType::k01s);



    py::class_<SlipBondType> slip_bond_t(m, "SlipBondType", abstr_bond_t);

    slip_bond_t.def(
        py::init<double, double, double, double, double, double>(),
        "lambda_"_a, "sigma"_a, "k_f_0"_a, "rec_dens"_a, "x1s"_a, "k01s"_a
    );

    py::class_<CatchSlipBondType> cs_bond_type(m, "__CatchSlipBondType", abstr_bond_t);

    cs_bond_type.def_readwrite("x1c", &CatchSlipBondType::x1c)
                .def_readwrite("k01c", &CatchSlipBondType::k01c);

    py::class_<CatchSlipPselBondType> cs_psel_bond_t(m, "CatchSlipPselBondType", cs_bond_type);

    cs_psel_bond_t.def(
            py::init<double, double, double, double, double, double, double, double>(),
            "lambda_"_a, "sigma"_a, "k_f_0"_a, "rec_dens"_a, "x1s"_a, "k01s"_a, "x1c"_a, "k01c"_a
    );

    py::class_<CatchSlipIntegrinBondType> cs_int_bond_t(m, "CatchSlipIntegrinBondType", cs_bond_type);

    cs_int_bond_t.def(
            py::init<double, double, double, double, double, double, double, double>(),
            "lambda_"_a, "sigma"_a, "k_f_0"_a, "rec_dens"_a, "x1s"_a, "k01s"_a, "x1c"_a, "k01c"_a
    );


    ////////////////
    //  Settings  //
    ////////////////

    py::class_<Settings> settings(m, "Settings");

    settings.def(py::init<Settings::ModelParameters*>(), "p"_a)
            .def("add_lig_type", &Settings::add_lig_type, "lig_type"_a, "n_of_lig"_a);

    settings.def_readwrite("p", &Settings::p)
            .def_readwrite("lig_types_and_nrs", &Settings::lig_types_and_nrs);

    // ModelParameters

    py::class_<Settings::ModelParameters> p(settings, "ModelParameters");

    p.def(
            py::init<double, double, double, double, double, double>(),
            "r_c"_a, "mu"_a, "temp"_a, "dens_diff"_a, "f_rep_0"_a, "tau"_a
    );

    p.def_readwrite("r_c", &Settings::ModelParameters::r_c)
            .def_readwrite("mu", &Settings::ModelParameters::mu)
            .def_readwrite("temp", &Settings::ModelParameters::temp)
            .def_readwrite("dens_diff", &Settings::ModelParameters::dens_diff)
            .def_readwrite("f_rep_0", &Settings::ModelParameters::f_rep_0)
            .def_readwrite("tau", &Settings::ModelParameters::tau);

    // LigandType

    py::class_<Settings::LigandType> lig_type(settings, "LigandType");

    // passing a pointer is OK
    lig_type.def(py::init())
            .def("add_bond_type", &Settings::LigandType::add_bond_type, "bond_type"_a);

    lig_type.def_readonly("bonds_types", &Settings::LigandType::bonds_types);


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


