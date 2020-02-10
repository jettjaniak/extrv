#pragma once

#include "types.h"

namespace helpers {
    /**
     * Compute (r_cir, alpha) parametrization from (x, y) coordinates of point laying on a sphere.
     * @param lig_x x ligand coordinate
     * @param lig_y y ligand coordinate
     * @return (r_cir, alpha), where alpha is in [-π, π]
     */
    pair<double, double> parametrize_ligand(double lig_x, double lig_y);

    /**
     * Compute rate of binding reaction in Bell model.
     *
     * @param deviation absolute deviation from optimal bond length in μm
     * @param rate_0 reaction rate when deviation is zero in 1/s
     * @param spring_const spring constant in kg/s^2
     * @param react_compl reactive compliance in μm
     * @return bonding rate in 1/s
     */
    double bell_binding_rate(double deviation, double rate_0, double spring_const, double react_compl);

    /**
     * Draw from Uniform([0,1]) distribution.
     *
     * @param generator random number generator
     * @return number in [0, 1]
     */
    double draw_from_uniform_dist(generator_t generator);
}