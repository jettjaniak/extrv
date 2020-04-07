#pragma once
#include "AbstrConstStepSS.h"
#include <boost/numeric/odeint/stepper/euler.hpp>
using namespace boost::numeric::odeint;

typedef euler<array<double, 3>> euler_stepper_type;

struct AbstrEulerSS : AbstrConstStepSS {
    euler_stepper_type stepper;

    explicit AbstrEulerSS(double dt);

    void reset_stepper() override;
    double do_ode_step() override;
};

struct EulerGillSS : AbstrEulerSS, AbstrGillSS {
    EulerGillSS(double h_0, Parameters* p, unsigned int seed, double dt);
};

struct EulerProbSS : AbstrEulerSS, AbstrProbSS {
    EulerProbSS(double h_0, Parameters* p, unsigned int seed, double dt);
};