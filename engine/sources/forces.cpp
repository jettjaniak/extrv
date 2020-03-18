#include "forces.h"

#include "interpolated.h"


namespace forces {
    forces_t non_bond_forces(double shear_rate, double h, const Parameters *p) {
        forces_t f;
        f.f_y = nonspecific_force(h, p) + grav_force(p);
        if (shear_rate != 0.0)
            f += shear_forces(shear_rate, h, p);
        return f;
    }

    forces_t shear_forces(double shear_rate, double h, const Parameters *p) {
        double f_x = 6 * PI * p->visc * p->r_cell * (p->r_cell + h) * shear_rate * interpolated::f_s(h / p->r_cell);
        double t_z = -4 * PI * p->visc * shear_rate * interpolated::t_s(h / p->r_cell);
        return {f_x, 0.0, t_z};
    }

    double grav_force(const Parameters *p) {
        return -(4.0 / 3.0) * PI * pow(p->r_cell, 3.0) * p->dens_diff * G;
    }

    double van_der_waals_potential(double h, const Parameters *p) {
        static const double A_H = 5e-21 * 1e12;
        static const double CHI = 70 * 1e-4;
        
        double w_1 = 1 / pow(h, 2);
        double w_2 = -2 / pow(h + CHI, 2);
        double w_3 = 1 / pow(h + 2 * CHI, 2);
        return -(A_H / (12 * PI)) * (w_1 + w_2 + w_3);
    }

    double electrostatic_potential(double h, const Parameters *p) {
        static const double RHO = 5 * 1e-12;
        static const double KAPPA = (1.0 / 8) * 1e4;
        static const double L_C = 10 * 1e-3;
        static const double EPS = 7.08 * 1e-19 * 1e-9;
        
        double rho2 = pow(RHO, 2);
        double kappa3 = pow(KAPPA, 3);
        if (h > 2 * L_C) {
            // (A12)
            // in the 1984 reference we have kappa^2 in the denominator,
            // but it isn't consistent with the plots
            double f_e_1 = rho2 / (2 * EPS * kappa3);
            double f_e_2 = pow(exp(KAPPA * L_C) - 1, 2) * exp(-KAPPA * h);
            return f_e_1 * f_e_2;
        } else {
            // (A13)
            double f_e_1 = (rho2 / (2 * EPS * pow(KAPPA, 2))) * (4 * L_C - 2 * h);

            double f_e_2_1 = (rho2 / (2 * EPS * kappa3)) * exp(-KAPPA * h);
            double f_e_2_2 = (1 - 2 * exp(KAPPA * L_C));
            double f_e_2 = f_e_2_1 * f_e_2_2;

            // in the 1984 reference we have kappa^e in the denominator,
            // but it isn't consistent with the plots
            double f_e_3_1 = (rho2 / (2 * EPS * kappa3)) * exp(KAPPA * h);
            double f_e_3_2 = exp(-2 * KAPPA * L_C);
            double f_e_3 = f_e_3_1 * f_e_3_2;

            return f_e_1 + f_e_2 + f_e_3;
        }
    }

    double steric_stabilization_potential(double h, const Parameters *p) {
        static const double LAM_SS = 2.5e-6 * 10;
        static const double L_C = 10 * 1e-3;

        if (h > 2 * L_C)
            return 0;  // (A14)
        else
            return LAM_SS * ((2 * L_C - h) / pow(L_C, 2)); // (A15)
    }

    double nonspecific_force(double h, const Parameters *p) {
        double potential = van_der_waals_potential(h, p) + electrostatic_potential(h, p) +
                           steric_stabilization_potential(h, p);
        return 2 * PI * p->r_cell * potential;
    }
}