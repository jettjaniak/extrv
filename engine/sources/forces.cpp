#include "forces.h"

#include "interpolated.h"


forces_t forces::non_bond_forces(double shear_rate, double h, const Parameters *p) {
    forces_t f;
    f.f_y = repulsive_force(h, p) + grav_force(p);
    if (shear_rate != 0.0)
        f += shear_forces(shear_rate, h, p);
    return f;
}

forces_t forces::shear_forces(double shear_rate, double h, const Parameters *p) {
    double f_x = 6 * PI * p->visc * p->r_cell * (p->r_cell + h) * shear_rate * interpolated::f_s(h / p->r_cell);
    double t_z = -4 * PI * p->visc * shear_rate * interpolated::t_s(h / p->r_cell);
    return {f_x, 0.0, t_z};
}

double forces::repulsive_force(double h, const Parameters *p) {
    double f_rep_0 = p->f_rep_0;
    double tau = p->tau;
    // formula is constructed for Å, hence * 1e4
    double exp_m_tau_h = exp(-tau * h * 1e4);
    return f_rep_0 * (tau * exp_m_tau_h) / (1.0 - exp_m_tau_h);
}

double forces::grav_force(const Parameters *p) {
    return -(4.0 / 3.0) * PI * pow(p->r_cell, 3.0) * p->dens_diff * G;
}
