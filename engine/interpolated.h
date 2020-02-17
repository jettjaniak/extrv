#pragma once

#include "types.h"
#include "helpers.h"

#include <algorithm>


/**
 * Functions related to fluid dynamic, interpolated from data tables.
 *
 * All of them are similar to exponential function, so we linearly interpolate
 * logarithms of their values and then apply exp.
 *
 * Sometimes asymptotic formulas are used for small argument values.
 */
namespace interpolated {

    double f_r(double h_over_r);

    double f_s(double h_over_r);

    double f_t(double h_over_r);

    double t_r(double h_over_r);

    double t_s(double h_over_r);

    double t_t(double h_over_r);

    // ignore static storage initialization warnings, begin
    #pragma clang diagnostic push
    #pragma ide diagnostic ignored "cert-err58-cpp"

    const vector<double> F_R_DATA_X {3.2020e-03, 5.0040e-03, 4.5300e-02, 1.2760e-01,
                                     5.4310e-01, 1.3524e+00, 2.7622e+00, 9.0677e+00};

    const vector<double> F_R_DATA_Y {5.1326e-01, 4.5582e-01, 1.9403e-01, 9.8291e-02,
                                     1.9532e-02, 3.5231e-03, 5.6214e-04, 1.1699e-05};

    const vector<double> LOG_F_R_DATA_Y = helpers::transform_vector(F_R_DATA_Y, helpers::log_dbl);

    const vector<double> F_S_DATA_X {0.0000e+00, 3.2020e-03, 5.0040e-03, 4.5300e-02, 1.2760e-01,
                                     5.4310e-01, 1.3524e+00, 2.7622e+00, 9.0677e+00};

    const vector<double> F_S_DATA_Y {1.7005e+00, 1.6982e+00, 1.6969e+00, 1.6682e+00, 1.6160e+00,
                                     1.4391e+00, 1.2780e+00, 1.1671e+00, 1.0587e+00};

    const vector<double> LOG_F_S_DATA_Y_MINUS_1 = helpers::transform_vector(F_S_DATA_Y, helpers::log_x_minus_one);


    const vector<double> F_T_DATA_X {3.2020e-03, 5.0040e-03, 4.5300e-02, 1.2760e-01, 5.4310e-01,
                                     1.3524e+00, 2.7622e+00, 9.0677e+00};

    const vector<double> F_T_DATA_Y {-4.0223e+00, -3.7863e+00, -2.6475e+00, -2.1514e+00, -1.5675e+00,
                                     -1.3097e+00, -1.1738e+00, -1.0591e+00};

    const vector<double> LOG_MINUS_F_T_DATA_Y_MINUS_1 =
            helpers::transform_vector(F_T_DATA_Y, helpers::log_minus_x_minus_one);

    const vector<double> T_R_DATA_X {3.2020e-03, 5.0040e-03, 4.5300e-02, 1.2760e-01, 5.4310e-01,
                                     1.3524e+00, 2.7622e+00, 9.0677e+00};

    const vector<double> T_R_DATA_Y {-2.6793e+00, -2.5056e+00, -1.6996e+00, -1.3877e+00, -1.0998e+00,
                                     -1.0250e+00, -1.0059e+00, -1.0003e+00};

    const vector<double> LOG_MINUS_T_R_DATA_Y_MINUS_1 =
            helpers::transform_vector(T_R_DATA_Y, helpers::log_minus_x_minus_one);

    const vector<double> T_S_DATA_X {0.0000e+00, 3.2020e-03, 5.0040e-03, 4.5300e-02, 1.2760e-01,
                                     5.4310e-01, 1.3524e+00, 2.7622e+00, 9.0677e+00};

    const vector<double> T_S_DATA_Y {9.4399e-01, 9.4427e-01, 9.4442e-01, 9.4769e-01, 9.5374e-01,
                                     9.7419e-01, 9.9010e-01, 9.9711e-01, 9.9981e-01};

    const vector<double> LOG_MINUS_T_S_DATA_Y_PLUS_1 =
            helpers::transform_vector(T_S_DATA_Y, helpers::log_minus_x_plus_one);

    const vector<double> T_T_DATA_X {3.2020e-03, 5.0040e-03, 4.5300e-02, 1.2760e-01, 5.4310e-01,
                                     1.3524e+00, 2.7622e+00, 9.0677e+00};

    const vector<double> T_T_DATA_Y {3.8494e-01, 3.4187e-01, 1.4552e-01, 7.3718e-02, 1.4649e-02,
                                     2.6423e-03, 4.2160e-04, 8.7744e-06};

    const vector<double> LOG_T_T_DATA_Y = helpers::transform_vector(T_T_DATA_Y, helpers::log_dbl);

    // ignore static storage initialization warnings, end
    #pragma clang diagnostic pop
}

