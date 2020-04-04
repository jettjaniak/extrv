#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// for string formatting
#include <iomanip>
#include <sstream>

#include "types.h"
#include "Parameters.h"
#include "AbstractSimulationState.h"
#include "AbstractConstStepSimulationState.h"
#include "EulerSimulationState.h"
#include "RKSimulationState.h"
#include "AdaptiveSimulationState.h"
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
            py::init<double, double, double, double, double, double>(),
            "r_cell"_a, "visc"_a, "temp"_a, "dens_diff"_a,
            "rep_0"_a, "rep_scale"_a
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


    ///////////////////////////////
    //  AbstractSimulationState  //
    ///////////////////////////////

    py::class_<AbstractSimulationState> a_ss(m, "AbstractSimulationState");

    // Diagnostic
    py::class_<AbstractSimulationState::Diagnostic> diag(a_ss, "Diagnostic");

    diag.def_readonly("dt_freq", &AbstractSimulationState::Diagnostic::dt_freq)
        .def_readonly("n_pos_not_ok", &AbstractSimulationState::Diagnostic::n_pos_not_ok);

    a_ss.def("simulate_one_step", &AbstractSimulationState::simulate_one_step)
        .def("simulate", &AbstractSimulationState::simulate,
            "max_time"_a, "max_steps"_a)
        .def("simulate_with_history", &AbstractSimulationState::simulate_with_history,
            "max_time"_a, "max_steps"_a, "save_every"_a=1e-2);

    a_ss.def_readwrite("time", &AbstractSimulationState::time)
        .def_readwrite("pos", &AbstractSimulationState::pos)
        .def_readwrite("bd_lig_ind", &AbstractSimulationState::bd_lig_ind)
        .def_readwrite("p", &AbstractSimulationState::p)
        .def_readwrite("shear_rate", &AbstractSimulationState::shear_rate)
        .def_readwrite("try_dt", &AbstractSimulationState::try_dt)
        .def_readonly("diag", &AbstractSimulationState::diag)
        .def_readonly("n_active_lig", &AbstractSimulationState::n_active_lig);

    // AbstractConstStepSimulationState
    py::class_<AbstractConstStepSimulationState> acs_ss(m, "AbstractConstStepSimulationState", a_ss);
    acs_ss.def_readwrite("dt", &AbstractConstStepSimulationState::dt);

    // EulerSimulationState
    py::class_<EulerSimulationState> euler_ss(m, "EulerSimulationState", acs_ss);
    euler_ss.def(
        py::init<double, Parameters*, unsigned int, double>(),
        "h_0"_a, "p"_a, "seed"_a, "dt"_a
    );

    // RKSimulationState
    py::class_<RKSimulationState> rk_ss(m, "RKSimulationState", acs_ss);
    rk_ss.def(
            py::init<double, Parameters*, unsigned int, double>(),
            "h_0"_a, "p"_a, "seed"_a, "dt"_a
    );

    // AdaptiveSimulationState
    py::class_<AdaptiveSimulationState> adap_ss(m, "AdaptiveSimulationState", a_ss);
    adap_ss.def_readwrite("max_dt", &AdaptiveSimulationState::max_dt);
    adap_ss.def(
            py::init<double, Parameters*, unsigned int, double, double, double>(),
            "h_0"_a, "p"_a, "seed"_a,
            "max_dt"_a = 0.1, "abs_err"_a=1e-10, "rel_err"_a=1e-6
    );
}


