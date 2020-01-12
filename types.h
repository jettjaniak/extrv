#pragma once

#include <vector>
#include <set>
#include <random>

#include "constants.h"

using std::vector;
using std::set;
using std::pair;


struct forces_t {
    double f_x = 0.0;
    double f_y = 0.0;
    double t_z = 0.0;

    forces_t& operator+=(const forces_t& f2) {
        f_x += f2.f_x;
        f_y += f2.f_y;
        t_z += f2.t_z;
        return *this;
    }
};

struct velocities_t {
    double v_x;
    double v_y;
    double o_z;
};

typedef std::default_random_engine generator_t;