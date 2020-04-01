#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// for string formatting
#include <iomanip>
#include <sstream>

#include "types.h"
#include "Parameters.h"
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


    abstr_bond_t.def_readwrite("eq_bond_len", &AbstractBondType::eq_bond_len)
                .def_readwrite("spring_const", &AbstractBondType::spring_const)
                .def_readwrite("binding_rate_0", &AbstractBondType::binding_rate_0)
                .def_readwrite("rec_dens", &AbstractBondType::rec_dens)
                .def_readwrite("react_compl_slip", &AbstractBondType::react_compl_slip)
                .def_readwrite("rup_rate_0_slip", &AbstractBondType::rup_rate_0_slip);



    py::class_<SlipBondType> slip_bond_t(m, "SlipBondType", abstr_bond_t);

    slip_bond_t.def(
        py::init<double, double, double, double, double, double>(),
        "eq_bond_len"_a, "spring_const"_a, "binding_rate_0"_a, "rec_dens"_a, "react_compl_slip"_a, "rup_rate_0_slip"_a
    );

    py::class_<AbstractCatchSlipBondType> cs_bond_type(m, "AbstractCatchSlipBondType", abstr_bond_t);

    cs_bond_type.def_readwrite("react_compl_catch", &AbstractCatchSlipBondType::react_compl_catch)
                .def_readwrite("rup_rate_0_catch", &AbstractCatchSlipBondType::rup_rate_0_catch);

    py::class_<CatchSlipPselBondType> cs_psel_bond_t(m, "CatchSlipPselBondType", cs_bond_type);

    cs_psel_bond_t.def(
            py::init<double, double, double, double, double, double, double, double>(),
            "eq_bond_len"_a, "spring_const"_a, "binding_rate_0"_a, "rec_dens"_a,
            "react_compl_slip"_a, "rup_rate_0_slip"_a, "react_compl_catch"_a, "rup_rate_0_catch"_a
    );

    py::class_<CatchSlipIntegrinBondType> cs_int_bond_t(m, "CatchSlipIntegrinBondType", cs_bond_type);

    cs_int_bond_t.def(
            py::init<double, double, double, double, double, double, double, double>(),
            "eq_bond_len"_a, "spring_const"_a, "binding_rate_0"_a, "rec_dens"_a,
            "react_compl_slip"_a, "rup_rate_0_slip"_a, "react_compl_catch"_a, "rup_rate_0_catch"_a
    );


    //////////////////
    //  Parameters  //
    //////////////////

    py::class_<Parameters> p(m, "Parameters");

    // LigandType
    // TODO: no docstring in Python
    py::class_<Parameters::LigandType> lig_type(p, "LigandType");

    // passing a pointer is OK
    lig_type.def(py::init())
            .def("add_bond_type", &Parameters::LigandType::add_bond_type, "bond_type"_a);

    lig_type.def_readonly("bonds_types", &Parameters::LigandType::bonds_types);


    p.def(
            py::init<double, double, double, double, double, double, double, double>(),
            "r_cell"_a, "visc"_a, "temp"_a, "dens_diff"_a,
            "rep_0"_a, "rep_scale"_a, "abs_err"_a, "rel_err"_a
    );

    p.def("add_ligands", &Parameters::add_ligands, "lig_type"_a, "n_of_lig"_a);

    p.def_readwrite("r_cell", &Parameters::r_cell)
     .def_readwrite("visc", &Parameters::visc)
     .def_readwrite("temp", &Parameters::temp)
     .def_readwrite("dens_diff", &Parameters::dens_diff)
     .def_readwrite("rep_0", &Parameters::rep_0)
     .def_readwrite("rep_scale", &Parameters::rep_scale)
     .def_readwrite("lig_types_and_nrs", &Parameters::lig_types_and_nrs);


    ///////////////
    //  History  //
    ///////////////

    py::class_<History> hist(m, "History");

    py::class_<History::BondTrajectory> bond_traj(hist, "BondTrajectory");

    bond_traj.def_readonly("start_i", &History::BondTrajectory::start_i)
             .def_readonly("positions", &History::BondTrajectory::positions);

    hist.def_readonly("time", &History::time)
        .def_readonly("h", &History::h)
        .def_readonly("rot", &History::rot)
        .def_readonly("dist", &History::dist)
        .def_readonly("bond_trajectories", &History::bond_trajectories);


    ///////////////////////
    //  SimulationState  //
    ///////////////////////

    py::class_<SimulationState> s(m, "SimulationState");

    // Diagnostic
    py::class_<SimulationState::Diagnostic> diag(s, "Diagnostic");

    diag.def_readonly("dt_freq", &SimulationState::Diagnostic::dt_freq)
        .def_readonly("n_pos_not_ok", &SimulationState::Diagnostic::n_pos_not_ok);

    s.def(
        py::init<double, Parameters*, unsigned int>(),
        "h_0"_a, "p"_a, "seed"_a
    );

    s.def("simulate_one_step", &SimulationState::simulate_one_step)
     .def("simulate", &SimulationState::simulate, "max_time"_a, "max_steps"_a)
     .def("simulate_with_history", &SimulationState::simulate_with_history,
          "max_time"_a, "max_steps"_a, "save_every"_a=1e-2);

    s.def_readwrite("pos", &SimulationState::pos)
     .def_readwrite("bd_lig_ind", &SimulationState::bd_lig_ind)
     .def_readwrite("p", &SimulationState::p)
     .def_readwrite("shear_rate", &SimulationState::shear_rate)
     .def_readonly("diag", &SimulationState::diag)
     .def_readonly("n_active_lig", &SimulationState::n_active_lig);
}


