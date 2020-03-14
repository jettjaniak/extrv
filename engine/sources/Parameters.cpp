#include "Parameters.h"

#include "helpers.h"


vector<double> Parameters::LigandType::binding_rates(double surface_dist) {
    vector<double> ret(bonds_types.size());
    AbstractBondType* bond_type;
    for (int i = 0; i < bonds_types.size(); i++) {
        bond_type = bonds_types[i];
        ret[i] = bond_type->binding_rate(surface_dist, p->temp);
    }
    return ret;
}

void Parameters::LigandType::add_bond_type(AbstractBondType *bond_type) {
    bonds_types.push_back(bond_type);
    max_surf_dist = compute_max_surf_dist();
}

double Parameters::LigandType::compute_max_surf_dist() {
    return helpers::bisection(1e-6, 1.0,
            [&](double sd) -> double {
                return max_binding_rate(sd) - MIN_RATE;
            });
}

double Parameters::LigandType::max_binding_rate(double surf_dist) {
    vector<double> rates = binding_rates(surf_dist);
    return *std::max_element(rates.begin(), rates.end());
}

Parameters::Parameters(double r_cell, double visc, double temp, double dens_diff, double f_rep_0, double tau) :
        r_cell(r_cell),
        // conversion to kg/(μm s)
        visc(visc * 1e-7),
        temp(temp),
        // conversion to kg/μm^3
        dens_diff(dens_diff * 1e-15),
        // conversion to kg μm / s^2
        f_rep_0(f_rep_0 * 1e-6),
        tau(tau) {}

void Parameters::add_ligands(LigandType *lig_type, size_t n_of_lig) {
    lig_type->p = this;
    lig_types_and_nrs.emplace_back(lig_type, n_of_lig);
    max_surf_dist = std::max(max_surf_dist, lig_type->max_surf_dist);
}
