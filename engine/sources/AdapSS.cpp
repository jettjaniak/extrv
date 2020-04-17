#include "AdapSS.h"

#include <boost/numeric/odeint/stepper/generation.hpp>


AbstrAdapSS::AbstrAdapSS(double max_dt, double max_dt_with_bonds, double abs_err, double rel_err) :
        max_dt(max_dt),
        max_dt_with_bonds(max_dt_with_bonds),
        abs_err(abs_err),
        rel_err(rel_err)

{
    reset_stepper();
    try_dt = max_dt;
}


double AbstrAdapSS::do_ode_step() {
    if (bd_lig_ind.empty())
        try_dt = std::min(try_dt, max_dt);
    else
        try_dt = std::min(try_dt, max_dt_with_bonds);
    // TODO: define as class parameter or define operator()
    namespace pl = std::placeholders;
    auto rhs_system = std::bind(&AbstrSS::rhs, std::ref(*this), pl::_1 , pl::_2 , pl::_3);

    controlled_step_result result;
    double dt_inout = 0.0, time_inout = 0.0;
    array<double, 3> pos_inout {};

    bool step_done = false;
    while (!step_done) {
        dt_inout = try_dt;
        pos_inout = pos;
        time_inout = time;
        result = stepper.try_step(rhs_system, pos_inout, time_inout, dt_inout);
        if (result == fail) {
            // solver made dt_inout smaller
            try_dt = dt_inout;
            continue;
        }
        step_done = true;
    }

    double step_done_with_dt = try_dt;
    try_dt = dt_inout;  // possibly larger after successful step
    pos = pos_inout;
    time = time_inout;

    return step_done_with_dt;
}

void AbstrAdapSS::reset_stepper() {
    stepper = make_controlled<error_stepper_type>(abs_err, rel_err);
}

AdapGillSS::AdapGillSS(
        double h_0, Parameters* p, unsigned int seed,
        double max_dt, double max_dt_with_bonds,
        double abs_err, double rel_err) :

        AbstrSS(h_0, p, seed),
        AbstrAdapSS(max_dt, max_dt_with_bonds, abs_err, rel_err),
        AbstrGillSS() {}

AdapProbSS::AdapProbSS(
        double h_0, Parameters* p, unsigned int seed,
        double max_dt, double max_dt_with_bonds,
        double abs_err, double rel_err) :

        AbstrSS(h_0, p, seed),
        AbstrAdapSS(max_dt, max_dt_with_bonds, abs_err, rel_err),
        AbstrProbSS() {}