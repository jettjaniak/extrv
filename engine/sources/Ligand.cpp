#include "Ligand.h"

#include "helpers.h"
#include "AbstractBondType.h"


// Public methods

Ligand::Ligand(xy_t lig_xy, Parameters::LigandType *lig_type) : lig_type(lig_type) {
    pair<double, double> r_alpha_pair = helpers::parametrize_ligand(lig_xy);
    r_cir = r_alpha_pair.first;
    rot_inc = r_alpha_pair.second;
    binding_rates.resize(lig_type->bonds_types.size());
}

bool Ligand::prepare_binding(double h, double rot, double dt, generator_t &generator) {
    if (bond_state != 0)  // TODO: why not > 0?
        return false;

    lig_type->compute_binding_rates(surface_dist(h, rot), binding_rates);

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

bool Ligand::prepare_rupture(double h, double rot, double dt, generator_t &generator) {
    AbstractBondType* bond_type = get_curr_bond_type();  // will abort if not bonded
    double rupture_rate = bond_type->rupture_rate(bond_length(h, rot), lig_type->p->temp);
    double rupture_probability = 1.0 - exp(-rupture_rate * dt);
    return helpers::draw_from_uniform_dist(generator) < rupture_probability;
}

void Ligand::bond(double rot) {
    bond_state = prepared_bond_state;
    prepared_bond_state = -1;
    bd_rec_x = x_pos(rot);
}

void Ligand::rupture() {
    bond_state = 0;
    bd_rec_x = INFTY;
}

void Ligand::move_bd_rec(double x_dist) {
    if (bd_rec_x == INFTY) abort();
    bd_rec_x -= x_dist;
}

forces_t Ligand::bond_forces(double h, double rot) {
    AbstractBondType* bond_type = get_curr_bond_type();

    double lig_x = x_pos(rot);
    double lig_y = y_pos(rot);

    xy_t bond_vector = helpers::compute_bond_vector(surface_dist(h, rot), lig_x, bd_rec_x);

    double bond_len = bond_vector.length();
    double f_common = (bond_type->spring_const * (bond_len - bond_type->eq_bond_len)) / bond_len;

    double f_x = f_common * bond_vector.x;
    double f_y = f_common * bond_vector.y;
    double t_z = f_y * lig_x - f_x * lig_y;

    return {f_x, f_y, t_z};
}

// Private methods

double Ligand::x_pos(double rot) const {
    return r_cir * sin(rot + rot_inc);
}

double Ligand::y_pos(double rot) const {
    return - r_cir * cos(rot + rot_inc);
}

double Ligand::surface_dist(double h, double rot) {
    double surf_dist = h + lig_type->p->r_cell + y_pos(rot);
    if (surf_dist < 0) abort();
    return surf_dist;
}

double Ligand::bond_length(double h, double rot) {
    if (bond_state < 1) abort();

    xy_t bond_vector = helpers::compute_bond_vector(surface_dist(h, rot), x_pos(rot), bd_rec_x);

    return bond_vector.length();
}

AbstractBondType* Ligand::get_curr_bond_type() {
    if (bond_state < 1 || bond_state > lig_type->bonds_types.size())
        abort();
    return lig_type->bonds_types[bond_state - 1];
}






