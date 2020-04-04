#include "RKSimulationState.h"


RKSimulationState::RKSimulationState(double h_0, Parameters* p, unsigned int seed, double dt) :
        AbstractConstStepSimulationState(h_0, p, seed, dt)
{
    reset_stepper();
}

void RKSimulationState::reset_stepper() {
    stepper = rk_stepper_type();
}

// TODO: define one level higher
double RKSimulationState::do_ode_step() {
    namespace pl = std::placeholders;
    auto rhs_system = std::bind(&AbstractSimulationState::rhs, std::ref(*this), pl::_1 , pl::_2 , pl::_3);
    stepper.do_step(rhs_system, pos, time, try_dt);
    double step_done_with_dt = try_dt;
    time += step_done_with_dt;  // do_step takes time by value
    try_dt = dt;
    return step_done_with_dt;
}
