#include "AdaptiveSimulationState.h"

#include <boost/numeric/odeint/stepper/generation.hpp>


AdaptiveSimulationState::AdaptiveSimulationState(
        double h_0, Parameters *p, unsigned int seed, double max_dt) :

        AbstractSimulationState(h_0, p, seed),
        max_dt(max_dt)

{
    reset_stepper();
    try_dt = max_dt;
}


double AdaptiveSimulationState::do_ode_step() {
    try_dt = std::min(try_dt, max_dt);
    // TODO: define as class parameter or define operator()
    namespace pl = std::placeholders;
    auto rhs_system = std::bind(&AbstractSimulationState::rhs, std::ref(*this), pl::_1 , pl::_2 , pl::_3);

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
    diag.add_dt(step_done_with_dt);  // diagnostics
    try_dt = dt_inout;  // possibly larger after successful step
    pos = pos_inout;
    time = time_inout;

    return step_done_with_dt;
}

void AdaptiveSimulationState::reset_stepper() {
    stepper = make_controlled<error_stepper_type>(abs_err, rel_err);
}
