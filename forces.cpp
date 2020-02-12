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
    double f_x = 6 * PI * p->mu * p->r_c * (p->r_c + h) * shear_rate * interpolated::f_s(h / p->r_c);
    double t_z = -4 * PI * p->mu * shear_rate * interpolated::t_s(h / p->r_c);
    return {f_x, 0.0, t_z};
}

double forces::repulsive_force(double h, const Parameters *p) {
    return 0; // TODO: implement
}

double forces::grav_force(const Parameters *p) {
    return 0; // TODO: implement
}
