#pragma once
#include "AbstractConstStepSimulationState.h"
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
using namespace boost::numeric::odeint;

typedef runge_kutta_dopri5<array<double, 3>> rk_stepper_type;


struct RKSimulationState : AbstractConstStepSimulationState {

    rk_stepper_type stepper;

    void reset_stepper() override;
    double do_ode_step() override;

    RKSimulationState(double h_0, Parameters* p, unsigned int seed, double dt);
};