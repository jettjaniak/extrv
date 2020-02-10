#include "Ligand.h"

#include "helpers.h"

#include <cmath>


Ligand::Ligand(double lig_x, double lig_y, LigandParameters* lig_p_) {
    auto r_alpha_pair = helpers::parametrize_ligand(lig_x, lig_y);
    r_cir = r_alpha_pair.first;
    alpha_inc = r_alpha_pair.second;

    lig_p = lig_p_;
}

double Ligand::surface_dist(double h, double alpha_0, Parameters *p) {
    return h + R_C + y_pos(alpha_0);
}

bool Ligand::prepare_binding(double h, double alpha_0, double dt, Parameters *p, generator_t generator) {
    if (bond_state != 0)
        return false;

    // TODO: more binding rates and bond states (different receptors)
    double deviation = std::abs(surface_dist(h, alpha_0, p) - lig_p->lambda_);
    double binding_rate = helpers::bell_binding_rate(deviation, lig_p->k_on_0, lig_p->sigma, lig_p->x1s);
    double binding_probability = 1.0 - exp(-binding_rate * dt);
    if (helpers::draw_from_uniform_dist(generator) < binding_probability) {
        prepared_bond_state = 1;
        return true;
    } else
        return false;

}

double Ligand::x_pos(double alpha_0) {
    return r_cir * sin(alpha_0 + alpha_inc);
}

double Ligand::y_pos(double alpha_0) {
    return r_cir * cos(alpha_0 + alpha_inc);
}

void Ligand::bond(double alpha_0) {
    bond_state = prepared_bond_state;
    prepared_bond_state = 0;
    bd_rec_x = x_pos(alpha_0);
}




