#include "../headers/Ligand.h"

#include "../headers/helpers.h"


// Public methods

Ligand::Ligand(xy_t lig_xy, LigandType *lig_type_, ModelParameters *p_) {
    pair<double, double> r_alpha_pair = helpers::parametrize_ligand(lig_xy);
    r_cir = r_alpha_pair.first;
    alpha_inc = r_alpha_pair.second;

    p = p_;
    lig_type = lig_type_;
}

bool Ligand::prepare_binding(double h, double alpha_0, double dt, generator_t &generator) {
    if (bond_state != 0)
        return false;

    vector<double> binding_rates = lig_type->binding_rates(surface_dist(h, alpha_0), p->temp);

    double any_binding_rate = 0.0;
    for (double binding_rate : binding_rates)
        any_binding_rate += binding_rate;

    double any_binding_probability = 1.0 - exp(-any_binding_rate * dt);
    if (helpers::draw_from_uniform_dist(generator) < any_binding_probability) {
        std::discrete_distribution<int> which_bond_distr {binding_rates.begin(), binding_rates.end()};
        prepared_bond_state = which_bond_distr(generator) + 1;
        return true;
    } else
        return false;

}

bool Ligand::prepare_rupture(double h, double alpha_0, double dt, generator_t &generator) {
    BondParameters* bond_p = get_curr_bond_p();  // will abort if not bonded
    double rupture_rate = bond_p->rupture_rate(bond_length(h, alpha_0), p->temp);
    double rupture_probability = 1.0 - exp(-rupture_rate * dt);
    return helpers::draw_from_uniform_dist(generator) < rupture_probability;
}

void Ligand::bond(double alpha_0) {
    bond_state = prepared_bond_state;
    prepared_bond_state = -1;
    bd_rec_x = x_pos(alpha_0);
}

void Ligand::rupture() {
    bond_state = 0;
    bd_rec_x = INFTY;
}

void Ligand::move_bd_rec(double x_dist) {
    if (bd_rec_x == INFTY) abort();
    bd_rec_x -= x_dist;
}

forces_t Ligand::bond_forces(double h, double alpha_0) {
    BondParameters* bond_p = get_curr_bond_p();

    double lig_x = x_pos(alpha_0);
    double lig_y = y_pos(alpha_0);

    xy_t bond_vector = helpers::compute_bond_vector(surface_dist(h, alpha_0), lig_x, bd_rec_x);

    double bond_len = bond_vector.length();
    double f_common = (bond_p->sigma * (bond_len - bond_p->lambda_)) / bond_len;

    double f_x = f_common * bond_vector.x;
    double f_y = f_common * bond_vector.y;
    double t_z = f_y * lig_x - f_x * lig_y;

    return {f_x, f_y, t_z};
}

// Private methods

double Ligand::x_pos(double alpha_0) {
    return r_cir * sin(alpha_0 + alpha_inc);
}

double Ligand::y_pos(double alpha_0) {
    return - r_cir * cos(alpha_0 + alpha_inc);
}

double Ligand::surface_dist(double h, double alpha_0) {
    double surf_dist = h + p->r_c + y_pos(alpha_0);
    if (surf_dist < 0) abort();
    return surf_dist;
}

double Ligand::bond_length(double h, double alpha_0) {
    if (bond_state < 1) abort();

    xy_t bond_vector = helpers::compute_bond_vector(surface_dist(h, alpha_0), x_pos(alpha_0), bd_rec_x);

    return bond_vector.length();
}

BondParameters* Ligand::get_curr_bond_p() {
    if (bond_state < 1 || bond_state > lig_type->bonds_p.size())
        abort();
    return lig_type->bonds_p[bond_state - 1];
}






