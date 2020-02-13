#include "helpers.h"


namespace helpers {
    pair<double, double> parametrize_ligand(xy_t lig_xy) {
        double r_cir = lig_xy.length();
        double alpha;
        // x = r sin(alpha)
        if (lig_xy.y < 0)
            // We are in the bottom half, where alpha is in [-π/2, π/2],
            // which corresponds to asin return value range.
            alpha = asin(lig_xy.x / r_cir);
        else {
            // We are in the top half, where alpha is in
            // [-π, -π/2] or [π/2, π].

            // sin(alpha - π) = - x / r
            alpha = PI - asin(lig_xy.x / r_cir);
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
        return rate_0 * exp(react_compl * force);
    }

    double psel_rupture_rate(double force, double rate_0_slip, double rate_0_catch, double react_compl_slip,
                             double react_compl_catch) {
        double half_f = force / 2;
        return (2.0 / 3.0) * rate_0_slip * exp(react_compl_slip * half_f)
            + rate_0_catch * exp(react_compl_slip * half_f);
    }

    double integrin_rupture_rate(double force, double rate_0_slip, double rate_0_catch, double react_compl_slip,
                                 double react_compl_catch) {
        return 1 / (
                (1 / rate_0_slip) * std::exp(-react_compl_slip * force) +
                (1 / rate_0_catch) * std::exp(-react_compl_catch * force)
        );
    }

    std::uniform_real_distribution<double> uniform_dist {};

    double draw_from_uniform_dist(generator_t &generator) {
        return uniform_dist(generator);
    }

    std::normal_distribution<double> normal_dist {};

    double draw_from_normal_dist(generator_t &generator) {
        return normal_dist(generator);
    }

    xyz_t draw_from_uniform_dist_on_sphere(double radius, generator_t &generator) {
        xyz_t ret {
            draw_from_normal_dist(generator),
            draw_from_normal_dist(generator),
            draw_from_normal_dist(generator)
        };

        return ret * (radius / ret.length());
    }

    xy_t compute_bond_vector(double surf_dist, double lig_x, double rec_x) {
        return {rec_x - lig_x, - surf_dist};
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

    double log_x_minus_one(double x) {
        return log(x - 1);
    }

    double log_minus_x_minus_one(double x) {
        return log(- x - 1);
    }

    double log_minus_x_plus_one(double x) {
        return log(- x + 1);
    }

}