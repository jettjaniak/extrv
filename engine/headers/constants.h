#pragma once

#include <limits>


constexpr double PI = 3.141592653589793238463;
/// Boltzmann constant in (μm^2 kg) / (s^2 K)
constexpr double K_B = 1.38064852e-11;
/// gravitational acceleration in μm / s^2 (in 1992 paper it's called gravitational constant)
constexpr double G = 9806650.;
constexpr double INFTY = std::numeric_limits<double>::infinity();
constexpr int LAMBDA_SERIES_MAX_N = 50;
constexpr double MIN_RATE = 1e-10;

constexpr int POS_LOG_H = 0;
constexpr int POS_ROT = 1;
constexpr int POS_DIST = 2;

constexpr double MAX_DT = 0.1;
constexpr double MAX_LOG_H = 2.5;

constexpr double DOUBLE_MIN = std::numeric_limits<double>::min();
constexpr double DOUBLE_DENORM_MIN = std::numeric_limits<double>::denorm_min();


