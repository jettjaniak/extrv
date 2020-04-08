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


void Ligand::bond(double rot, double dist, generator_t &generator) {
    std::discrete_distribution<int>
        which_bond_distr {binding_rates.begin(), binding_rates.end()};
    bond_state = which_bond_distr(generator) + 1;
    bd_rec_x = x_pos(rot) + dist;
}

void Ligand::rupture() {
    bond_state = 0;
    bd_rec_x = INFTY;
}

forces_t Ligand::bond_forces(const array<double, 3> & pos) {
    AbstractBondType* bond_type = get_curr_bond_type();

    double lig_x = x_pos(pos[POS_ROT]);
    double lig_y = y_pos(pos[POS_ROT]);

    xy_t bond_vector = helpers::compute_bond_vector(
            surface_dist(pos),
            lig_x, bd_rec_x - pos[POS_DIST]
    );

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

double Ligand::surface_dist(const array<double, 3> & pos) {
    double surf_dist = exp(pos[POS_LOG_H]) + lig_type->p->r_cell + y_pos(pos[POS_ROT]);
    if (surf_dist < 0)
        abort();
    return surf_dist;
}

double Ligand::bond_length(const array<double, 3> & pos) {
    if (bond_state < 1)
        abort();

    xy_t bond_vector = helpers::compute_bond_vector(
            surface_dist(pos), x_pos(pos[POS_ROT]), bd_rec_x - pos[POS_DIST]);

    return bond_vector.length();
}

AbstractBondType* Ligand::get_curr_bond_type() {
    if (bond_state < 1 || bond_state > lig_type->bonds_types.size())
        abort();
    return lig_type->bonds_types[bond_state - 1];
}

double Ligand::update_binding_rates(const array<double, 3> & pos) {
    double any_binding_rate = 0.0;
    AbstractBondType* bond_type;
    for (int i = 0; i < lig_type->bonds_types.size(); i++) {
        bond_type = lig_type->bonds_types[i];
        binding_rates[i] = bond_type->binding_rate(surface_dist(pos), lig_type->p->temp);
        any_binding_rate += binding_rates[i];
    }
    return any_binding_rate;
}

double Ligand::rupture_rate(const array<double, 3> & pos) {
    AbstractBondType* bond_type = get_curr_bond_type();  // will abort if not bonded
    return bond_type->rupture_rate(bond_length(pos), lig_type->p->temp);
}






