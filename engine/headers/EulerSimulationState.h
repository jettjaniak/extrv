#pragma once
#include "AbstractConstStepSimulationState.h"
#include <boost/numeric/odeint/stepper/euler.hpp>
using namespace boost::numeric::odeint;

typedef euler<array<double, 3>> euler_stepper_type;


struct EulerSimulationState : AbstractConstStepSimulationState {

    euler_stepper_type stepper;

    void reset_stepper() override;
    double do_ode_step() override;

    EulerSimulationState(double h_0, Parameters* p, unsigned int seed, double dt);
};