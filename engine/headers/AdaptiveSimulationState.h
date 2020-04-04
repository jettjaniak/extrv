#pragma once

#include "AbstractSimulationState.h"

#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>

using namespace boost::numeric::odeint;


typedef runge_kutta_dopri5<array<double, 3>> error_stepper_type;
typedef controlled_runge_kutta<error_stepper_type> adaptive_stepper_type;


struct AdaptiveSimulationState : AbstractSimulationState {

    adaptive_stepper_type stepper;
    double abs_err = 1e-10;
    double rel_err = 1e-6;
    double max_dt;

    void reset_stepper() override;
    double do_ode_step() override;

    AdaptiveSimulationState(double h_0, Parameters* p, unsigned int seed, double max_dt = 0.1);
};