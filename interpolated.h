#pragma once

#include "types.h"

/**
 * Functions related to fluid dynamic, interpolated from data tables.
 *
 * All of them are similar to exponential function, so we lineary interpolate
 * logarithms of their values and then apply exp.
 *
 * Sometimes asymptotic formulas are used for small argument values.
 */
// TODO: implement
namespace interpolated {

    double f_r(double h_over_r);

    double f_s(double h_over_r);

    double f_t(double h_over_r);

    double t_r(double h_over_r);

    double t_s(double h_over_r);

    double t_t(double h_over_r);

}