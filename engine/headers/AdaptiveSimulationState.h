#pragma once

#include "AbstractSimulationState.h"

#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>

using namespace boost::numeric::odeint;


typedef runge_kutta_dopri5<array<double, 3>> error_stepper_type;
typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;


struct AdaptiveSimulationState : AbstractSimulationState {

    controlled_stepper_type stepper;
    double abs_err = 1e-10;
    double rel_err = 1e-6;

    void reset_stepper() override;
    double do_ode_step() override;

    AdaptiveSimulationState(double h_0, Parameters* p, unsigned int seed);
};