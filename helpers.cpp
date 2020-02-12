#include "helpers.h"

#include <cmath>


namespace helpers {
    pair<double, double> parametrize_ligand(double lig_x, double lig_y) {
        double r_cir = sqrt(lig_x * lig_x + lig_y * lig_y);
        double alpha;
        // x = r sin(alpha)
        if (lig_y < 0)
            // We are in the bottom half, where alpha is in [-π/2, π/2],
            // which corresponds to asin return value range.
            alpha = asin(lig_x / r_cir);
        else {
            // We are in the top half, where alpha is in
            // [-π, -π/2] or [π/2, π].

            // sin(alpha - π) = - x / r
            alpha = PI - asin(lig_x / r_cir);
            // but we want alpha to stay in [-π, π]
            if (alpha > PI)
                alpha -= 2.0 * PI;
        }
        return {r_cir, alpha};
    }

    double bell_binding_rate(double deviation, double rate_0, double spring_const, double react_compl, double temp) {
        return rate_0 * exp((spring_const * deviation * (react_compl - 0.5 * deviation)) / (K_B * temp));
    }

    double esel_rupture_rate(double force, double rate_0, double react_compl) {
        return 0;  // TODO: implement
    }

    double psel_rupture_rate(double force, double rate_0_slip, double rate_0_catch, double react_compl_slip,
                             double react_compl_catch) {
        return 0;  // TODO: implement
    }

    double integrin_rupture_rate(double force, double rate_0_slip, double rate_0_catch, double react_compl_slip,
                                 double react_compl_catch) {
        return 0;  // TODO: implement
    }

    std::uniform_real_distribution<double> uniform_dist(0, 1);

    double draw_from_uniform_dist(generator_t generator) {
        return uniform_dist(generator);
    }

    pair<double, double> compute_bond_vector(double surf_dist, double lig_x, double rec_x) {
        return pair<double, double>();  // TODO: implement
    }

    double compute_2d_vector_length(pair<double, double> vector) {
        return 0;  // TODO: implement
    }

    /**
     * Right now implemented in dumb way with linear complexity.
     * Could be O(log n).
     */
    double linear_interpolation(const vector<double> &points_x, const vector<double> &points_y, double x) {
        if (x < points_x.front() || x > points_x.back())
            abort();

        size_t i = 1;
        while (x > points_x[i])
            i++;

        double y_0 = points_y[i - 1];
        double dist = x - points_x[i - 1];
        double slope = (points_y[i] - points_y[i - 1]) / (points_x[i] - points_x[i - 1]);

        return y_0 + dist * slope;
    }

}