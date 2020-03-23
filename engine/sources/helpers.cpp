#include "helpers.h"


namespace helpers {
    pair<double, double> parametrize_ligand(xy_t lig_xy) {
        double r_cir = lig_xy.length();  // by definition
        double rot_inc;
        // x = r sin(rot_inc)
        if (lig_xy.y < 0) {
            // We are in the bottom half, where rot_inc is in [3/2 π, 2π] or [0, π/2].
            // asin, which returns values in [-π/2, π/2] gives correct angle up to 2π
            rot_inc = asin(lig_xy.x / r_cir);                  // in [-π/2, π/2]
            rot_inc = std::fmod(rot_inc + 2 * PI, 2 * PI);  // in [0, 2π]
        } else {
            // We are in the top half, where rot_inc is in [π/2, 3/2 π]
            // sin(rot_inc - π) = - x / r
            rot_inc = PI - asin(lig_xy.x / r_cir);  // in [0, 2π]
        }
        return {r_cir, rot_inc};
    }

    double bell_binding_rate(double deviation, double rate_0, double spring_const, double react_compl, double temp) {
        return rate_0 * exp((spring_const * deviation * (react_compl - 0.5 * deviation)) / (K_B * temp));
    }

    double draw_from_uniform_dist(generator_t &generator) {
        static std::uniform_real_distribution<double> uniform_dist {};
        return uniform_dist(generator);
    }

    xyz_t draw_from_uniform_dist_on_sphere(double radius, generator_t &generator) {
        static std::normal_distribution<double> normal_dist {};
        xyz_t ret {
            normal_dist(generator),
            normal_dist(generator),
            normal_dist(generator)
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

    double log_dbl(double x) {
        return log(x);
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

    vector<double> transform_vector(const vector<double> &c, const std::function<double(double)> &f)
    {
        vector<double> ret;
        std::transform(std::begin(c), std::end(c), std::inserter(ret, std::end(ret)), f);
        return ret;
    }

    double bisection(double a, double b, const std::function<double(double)> &f, double tol) {
        double f_a = f(a);
        double f_b = f(b);
        if (f_a * f_b > 0 || a > b)
            abort();

        while (std::abs(f_b) > tol) {
            double new_point = a + (b - a) / 2;
            if (f(new_point) * f_b > 0) {
                b = new_point;
                f_b = f(b);
            }
            else
                a = new_point;
        }

        return b;
    }

}