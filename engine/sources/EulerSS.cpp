#include "EulerSS.h"


AbstrEulerSS::AbstrEulerSS(double dt) :
    AbstrConstStepSS(dt) {}

void AbstrEulerSS::reset_stepper() {}

double AbstrEulerSS::do_ode_step() {
    static array<double, 3> tiny_pos, dxdt;

//    static array<double, 3> dxdt;
    rhs(pos, dxdt, time);
//
//    for (int i = 0; i < 3; i++) {
//        pos[i] += dxdt[i] * try_dt;
//    }


    static double pos_item_temp, inc;
//    static long not_good = 0;
    for (int i = 0; i < 3; i++) {
        inc = tiny_pos[i] + dxdt[i] * try_dt;
        if (inc != 0.0) {
//            pos_item_temp = pos[i] + inc;

//            if (pos_item_temp != pos[i]) {
            // Good.
            if (pos[i] == 0.0 || std::abs(inc / pos[i]) > 1e-15) {
//                pos[i] = pos_item_temp;
                pos[i] += inc;
                tiny_pos[i] = 0.0;
            // Not good.
            } else {
                tiny_pos[i] = inc;
                // std::cout << ++not_good << std::endl;
            }
        }

//        else {
//            std::cout << "inc == 0.0";
//        }
    }


    namespace pl = std::placeholders;
    auto rhs_system = std::bind(&AbstrSS::rhs, std::ref(*this), pl::_1 , pl::_2 , pl::_3);

//    static array<double, 3> pos_out, dxdt;
//    stepper.do_step(rhs_system, pos, time, pos_out, try_dt);

//    vector<int> eq;
//    for (int i = 0; i < 3; i++) {
//        if (pos[i] == pos_out[i])
//            eq.push_back(i);
//    }
//    static long n_ebdn = 0;
//    if (!eq.empty()) {
//        static bool eq_but_dxdt_nonzero = false;
//        rhs(pos, dxdt, 0);
//        for (const auto &eq_ind : eq) {
//            if (dxdt[eq_ind] != 0) {
//                eq_but_dxdt_nonzero = true;
//                break;
//            }
//        }
//
//        if (eq_but_dxdt_nonzero) {
//            n_ebdn++;
//            std::cout << "< " << n_ebdn << "> " << std::endl;
//            std::cout << "pos before: [";
//            for (const auto &pos_item : pos) {
//                std::cout << pos_item << ", ";
//            }
//            std::cout << "]" << std::endl;
//
//            std::cout << "pos after:  [";
//            for (const auto &pos_item : pos_out) {
//                std::cout << pos_item << ", ";
//            }
//            std::cout << "]" << std::endl;
//
//            std::cout << "dx / dt:    [";
//            for (const auto &item : dxdt) {
//                std::cout << item << ", ";
//            }
//            std::cout << "]" << std::endl << std::endl;
//        }
//    }
//    pos = pos_out;


    double step_done_with_dt = try_dt;
    time += step_done_with_dt;  // do_step takes time by value
    try_dt = dt;
    return step_done_with_dt;
}

EulerGillSS::EulerGillSS(double h_0, Parameters* p, unsigned int seed, double dt) :
        AbstrSS(h_0, p, seed),
        AbstrEulerSS(dt),
        AbstrGillSS() {}

EulerProbSS::EulerProbSS(double h_0, Parameters* p, unsigned int seed, double dt) :
        AbstrSS(h_0, p, seed),
        AbstrEulerSS(dt),
        AbstrProbSS() {}