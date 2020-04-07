#include "RKSS.h"

AbstrRKSS::AbstrRKSS(double dt) : AbstrConstStepSS(dt) {
    reset_stepper();
}

void AbstrRKSS::reset_stepper() {
    stepper = rk_stepper_type();
}

double AbstrRKSS::do_ode_step() {
    namespace pl = std::placeholders;
    auto rhs_system = std::bind(&AbstrSS::rhs, std::ref(*this), pl::_1 , pl::_2 , pl::_3);
    stepper.do_step(rhs_system, pos, time, try_dt);
    double step_done_with_dt = try_dt;
    time += step_done_with_dt;  // do_step takes time by value
    try_dt = dt;
    return step_done_with_dt;
}


RKGillSS::RKGillSS(
        double h_0, Parameters *p, unsigned int seed, double dt) :

        AbstrSS(h_0, p, seed),
        AbstrRKSS(dt),
        AbstrGillSS() {}


RKProbSS::RKProbSS(
        double h_0, Parameters *p, unsigned int seed, double dt) :

        AbstrSS(h_0, p, seed),
        AbstrRKSS(dt),
        AbstrProbSS() {}
