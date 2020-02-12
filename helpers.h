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
     * @param temp temperature in K
     * @return binding rate in 1/s
     */
    double bell_binding_rate(double deviation, double rate_0, double spring_const, double react_compl, double temp);

    /**
     * Compute rate of rupture in slip E-selectin bond.
     *
     * @param force force exerted on bond in kg μm / s^2
     * @param rate_0 multiplicative constant in rupture rate in 1/s
     * @param react_compl reactive compliance in μm
     * @return rupture rate in 1/s
     */
    double esel_rupture_rate(double force, double rate_0, double react_compl);

    /**
     * Compute rate of rupture in catch-slip P-selectin bond.
     *
     * @param force force exerted on bond in kg μm / s^2
     * @param rate_0_slip multiplicative constant in slip part of rupture rate in 1/s
     * @param rate_0_catch multiplicative constant in catch part of rupture rate in 1/s
     * @param react_compl_slip reactive compliance in slip part of rupture rate in μm
     * @param react_compl_catch reactive compliance in catch part of rupture rate in μm
     * @return rupture rate in 1/s
     */
    double psel_rupture_rate(double force, double rate_0_slip, double rate_0_catch,
            double react_compl_slip, double react_compl_catch);

    /**
     * Compute rate of rupture in catch-slip ICAM and VCAM bonds.
     *
     * @param force force exerted on bond in kg μm / s^2
     * @param rate_0_slip multiplicative constant in slip part of rupture rate in 1/s
     * @param rate_0_catch multiplicative constant in catch part of rupture rate in 1/s
     * @param react_compl_slip reactive compliance in slip part of rupture rate in μm
     * @param react_compl_catch reactive compliance in catch part of rupture rate in μm
     * @return rupture rate in 1/s
     */
    double integrin_rupture_rate(double force, double rate_0_slip, double rate_0_catch,
            double react_compl_slip, double react_compl_catch);


    /**
     * Draw from Uniform([0,1]) distribution.
     *
     * @param generator random number generator
     * @return number in [0, 1]
     */
    double draw_from_uniform_dist(generator_t generator);

    /**
     * Compute vector of bond from ligand to corresponding receptor.
     *
     * @param surface_dist vertical distance from ligand to surface
	 * @param lig_x x coordinate of bonded ligand
	 * @param rec_x x coordinate of corresponding bonded receptor
     * @return
     */
    pair<double, double> compute_bond_vector(double surf_dist, double lig_x, double rec_x);

    /**
     * Compute Euclidian length.
     * @param vector x and y coordinates
     */
    double compute_2d_vector_length(pair<double, double> vector);
}