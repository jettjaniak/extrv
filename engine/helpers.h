#pragma once

#include "types.h"

#include <algorithm>

namespace helpers {
    /**
     * Compute (r_cir, alpha) parametrization from (x, y) coordinates of point laying on a sphere.
     * @param lig_xy x ligand coordinate
     * @param lig_y y ligand coordinate
     * @return (r_cir, alpha), where alpha is in [-π, π]
     */
    pair<double, double> parametrize_ligand(xy_t lig_xy);

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
    double esel_rupture_rate(double force, double rate_0, double react_compl, double temp);

    /**
     * Compute rate of rupture in catch-slip P-selectin bond.
     * TODO: citation
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
     * TODO: citation
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
     * Draw from Uniform([0, 1]) distribution.
     *
     * @param generator random number generator
     * @return number in [0, 1]
     */
    double draw_from_uniform_dist(generator_t &generator);

    /**
     * Draw from Normal(0, 1) distribution.
     *
     * @param generator random number generator
     * @return real number
     */
    double draw_from_normal_dist(generator_t &generator);

    /**
     * Draw from Uniform(S^2(0, radius)) distribution.
     *
     * @param generator random number generator
     * @param radius sphere radius
     * @return xyz triple
     */
    xyz_t draw_from_uniform_dist_on_sphere(double radius, generator_t &generator);

    /**
     * Compute vector of bond from ligand to corresponding receptor.
     *
     * @param surface_dist vertical distance from ligand to surface
	 * @param lig_x x coordinate of bonded ligand
	 * @param rec_x x coordinate of corresponding bonded receptor
     * @return
     */
    xy_t compute_bond_vector(double surf_dist, double lig_x, double rec_x);

    /**
     * Performs piece-wise linear interpolation.
     *
     * @param points arguments (ascending) in first column, values in second
     * @param x argument
     * @returns f(x), where f is piece-wise linear interpolant of points
     *
     */
    double linear_interpolation(const vector<double>& points_x, const vector<double>& points_y, double x);

    double log_x_minus_one(double x);
    double log_minus_x_minus_one(double x);
    double log_minus_x_plus_one(double x);



    /**
     * Returns container with elements transformed by functor.
     */
    template <typename Container, typename Functor>
    Container transform_container(const Container& c, Functor &&f)
    {
        Container ret;
        std::transform(std::begin(c), std::end(c), std::inserter(ret, std::end(ret)), f);
        return ret;
    }


}