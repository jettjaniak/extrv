#include "AdapSS.h"

#include <boost/numeric/odeint/stepper/generation.hpp>


AbstrAdapSS::AbstrAdapSS(double max_dt, double abs_err, double rel_err) :
        max_dt(max_dt),
        abs_err(abs_err),
        rel_err(rel_err)

{
    reset_stepper();
    try_dt = max_dt;
}


double AbstrAdapSS::do_ode_step() {
    try_dt = std::min(try_dt, max_dt);
    // TODO: define as class parameter or define operator()
    namespace pl = std::placeholders;
    auto rhs_system = std::bind(&AbstrSS::rhs, std::ref(*this), pl::_1 , pl::_2 , pl::_3);

    controlled_step_result result;
    double dt_inout = 0.0, time_inout = 0.0;
    array<double, 3> pos_inout {};

    bool dt_zero = false;
    bool step_done = false;
    while (!step_done) {
        dt_inout = try_dt;
        pos_inout = pos;
        time_inout = time;
        result = stepper.try_step(rhs_system, pos_inout, time_inout, dt_inout);
        // If solver is in the world of NaNs and infinities
        if (helpers::pos_not_ok(pos_inout)) {
            diag.n_pos_not_ok++;
            // we have to reset it
            reset_stepper();
            // and reduce the step size
            try_dt /= 2;
            continue;
        }
        if (result == fail) {
            // solver made dt_inout smaller
            try_dt = dt_inout;
            continue;
        }
        if (try_dt == 0) {
            if (dt_zero) {
                std::cout << "dt was 0 twice in a row, aborting" << std::endl;
                abort();
            } else {
                dt_zero = true;
            }
            try_dt = DOUBLE_MIN; // DOUBLE_DENORM_MIN;
            continue;
        }
        step_done = true;
    }

    double step_done_with_dt = try_dt;
    try_dt = std::min(dt_inout, max_dt);  // possibly larger after successful step
    pos = pos_inout;
    time = time_inout;

    return step_done_with_dt;
}

void AbstrAdapSS::reset_stepper() {
    stepper = make_controlled<error_stepper_type>(abs_err, rel_err);
}

AdapGillSS::AdapGillSS(
        double h_0, Parameters* p, unsigned int seed,
        double max_dt, double abs_err, double rel_err) :

        AbstrSS(h_0, p, seed),
        AbstrAdapSS(max_dt, abs_err, rel_err),
        AbstrGillSS() {}

AdapProbSS::AdapProbSS(
        double h_0, Parameters* p, unsigned int seed,
        double max_dt, double abs_err, double rel_err) :

        AbstrSS(h_0, p, seed),
        AbstrAdapSS(max_dt, abs_err, rel_err),
        AbstrProbSS() {}