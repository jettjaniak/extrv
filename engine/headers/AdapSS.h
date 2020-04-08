#pragma once

#include "AbstrSS.h"
#include "AbstrGillSS.h"
#include "AbstrProbSS.h"

#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>

using namespace boost::numeric::odeint;


typedef runge_kutta_dopri5<array<double, 3>> error_stepper_type;
typedef controlled_runge_kutta<error_stepper_type> adaptive_stepper_type;


struct AbstrAdapSS : virtual AbstrSS {

    adaptive_stepper_type stepper;

    /// maximal time step
    double max_dt;
    /// maximal time step
    double max_dt_with_bonds;
    /// limit on absolute error of ODE solver
    double abs_err;
    /// limit on relative error of ODE solver
    double rel_err;

    explicit AbstrAdapSS(double max_dt = 0.1, double max_dt_with_bonds = 1e-4,
            double abs_err = 1e-10, double rel_err = 1e-6);

    void reset_stepper() override;
    double do_ode_step() override;
};

struct AdapGillSS : AbstrAdapSS, AbstrGillSS {
    AdapGillSS(double h_0, Parameters* p, unsigned int seed,
               double max_dt = 0.1, double max_dt_with_bonds = 1e-4,
               double abs_err = 1e-10, double rel_err = 1e-6);
};

struct AdapProbSS : AbstrAdapSS, AbstrProbSS {
    AdapProbSS(double h_0, Parameters* p, unsigned int seed,
               double max_dt = 0.1, double max_dt_with_bonds = 1e-4,
               double abs_err = 1e-10, double rel_err = 1e-6);
};