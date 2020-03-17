#pragma once

#include "types.h"

#include <algorithm>
#include <functional>
#include <iterator>


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
     * Draw from Uniform([0, 1]) distribution.
     *
     * @param generator random number generator
     * @return number in [0, 1]
     */
    double draw_from_uniform_dist(generator_t &generator);

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

    double log_dbl(double x);  // there are some problems with std::function and overloaded functions in Visual C++
    double log_x_minus_one(double x);
    double log_minus_x_minus_one(double x);
    double log_minus_x_plus_one(double x);

    /**
     * Returns vector with elements transformed by f.
     */
    vector<double> transform_vector(const vector<double> &c, const std::function<double(double)> &f);

    /**
     * Find root of continuous function
     *
     * @param a left side of interval
     * @param b right side of interval
     * @param f function for which we want to find a root
     * @param tol tolerance
     * @return argument x of function for which |f(x)| < tol
     */
    double bisection(double a, double b, const std::function<double(double)> &f, double tol = 1e-6);
}