#pragma once
#include "AbstrConstStepSS.h"
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
using namespace boost::numeric::odeint;

typedef runge_kutta_dopri5<array<double, 3>> rk_stepper_type;

struct AbstrRKSS : AbstrConstStepSS {
    rk_stepper_type stepper;

    explicit AbstrRKSS(double dt);

    void reset_stepper() override;
    double do_ode_step() override;
};

struct RKGillSS : AbstrRKSS, AbstrGillSS {
    RKGillSS(double h_0, Parameters* p, unsigned int seed, double dt);
};

struct RKProbSS : AbstrRKSS, AbstrProbSS {
    RKProbSS(double h_0, Parameters* p, unsigned int seed, double dt);
};