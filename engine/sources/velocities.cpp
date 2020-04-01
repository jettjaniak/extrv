#include "velocities.h"

#include "interpolated.h"


namespace velocities {

    array<double, 3> compute_velocities(double log_h, const forces_t &f, const Parameters *p) {
        double h = exp(log_h);
        double h_over_r = h / p->r_cell;
        double t_t_val = interpolated::t_t(h_over_r);
        double f_r_val = interpolated::f_r(h_over_r);
        double f_t_val = interpolated::f_t(h_over_r);
        double t_r_val = interpolated::t_r(h_over_r);
        double d_24pi_mu_r2 = d_fun(t_t_val, f_r_val, f_t_val, t_r_val) * 24 * PI * p->visc * pow(p->r_cell, 2);
        double d_24pi_mu_r3 = d_24pi_mu_r2 * p->r_cell;
        double lambda = lambda_fun(h_over_r);

        // d dist / dt
        double v_x = (4 * p->r_cell * f.f_x * t_r_val + 3 * f_r_val * f.t_z) / d_24pi_mu_r2;
        // d h / dt
        double v_y = f.f_y / (6 * PI * lambda * p->visc * p->r_cell);
        // d rot / dt
        double o_z = (4 * p->r_cell * f.f_x * t_t_val + 3 * f_t_val * f.t_z) / d_24pi_mu_r3;

        array<double, 3> ret {};
        ret[POS_LOG_H] = v_y / h;  // d log(h) / dt = (1/h) * (d h / dt)
        ret[POS_ROT] = o_z;
        ret[POS_DIST] = v_x;

        return ret;
    }

    double lambda_fun(double h_over_r) {
        if (h_over_r < 1e-5)
            return 1 / h_over_r;

        double alpha = log(1.0 + h_over_r + sqrt(pow(1.0 + h_over_r, 2) - 1.0));
        double sum = 0.0;
        double p1, p2_num, p2_den, p2;
        // TODO: avoid large intermediate values
        for (int i = 1; i < LAMBDA_SERIES_MAX_N; i++) {
            p1 = (i * (i + 1.0)) / ((2.0 * i - 1) * (2.0 * i + 3.0));
            p2_num = 2.0 * sinh((2.0 * i + 1.0) * alpha) + (2.0 * i + 1.0) * sinh(2.0 * alpha);
            p2_den = 4.0 * pow(sinh((i + 0.5) * alpha), 2) - (2.0 * i + 1) * (2.0 * i + 1) * pow(sinh(alpha), 2);
            p2 = p2_num / p2_den - 1;
            sum += p1 * p2;
        }
        return (4.0 / 3.0) * sinh(alpha) * sum;
    }

    double d_fun(double t_t_val, double f_r_val, double f_t_val, double t_r_val) {
        return t_t_val * f_r_val - f_t_val * t_r_val;
    }

}