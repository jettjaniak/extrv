#pragma once

#include <limits>


constexpr double PI = 3.141592653589793238463;
/// Boltzmann constant in (μm^2 kg) / (s^2 K)
constexpr double K_B = 1.38064852e-11;
/// gravitational acceleration in μm / s^2 (in 1992 paper it's called gravitational constant)
constexpr double G = 9806650.;
constexpr double EPS_PROB = 1.1102230246251565e-16;
constexpr double INFTY = std::numeric_limits<double>::infinity();
constexpr int LAMBDA_SERIES_MAX_N = 50;