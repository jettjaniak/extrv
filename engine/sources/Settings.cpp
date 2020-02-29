#include "Settings.h"

#include "helpers.h"


Settings::ModelParameters::ModelParameters(double r_c, double mu, double temp, double dens_diff, double f_rep_0, double tau) :
        r_c(r_c),
        // conversion to kg/(μm s)
        mu(mu * 1e-7),
        temp(temp),
        // conversion to kg/μm^3
        dens_diff(dens_diff * 1e-15),
        // conversion to kg μm / s^2
        f_rep_0(f_rep_0 * 1e-6),
        tau(tau)
        {}

vector<double> Settings::LigandType::binding_rates(double surface_dist, double temp) {
    vector<double> ret(bonds_p.size());
    AbstractBondType* bond_p;
    for (int i = 0; i < bonds_p.size(); i++) {
        bond_p = bonds_p[i];
        ret[i] = bond_p->binding_rate(surface_dist, temp);
    }
    return ret;
}

void Settings::LigandType::add_bond_p(AbstractBondType *bond_p) {
    bonds_p.push_back(bond_p);
}

Settings::Settings(ModelParameters *p) : p(p) {}

void Settings::add_lig_type(LigandType *lig_type, size_t n_of_lig) {
    if (lig_type->index_in_settings > -1)
        abort();  // this ligand type was added to Settings previously, not implemented
    lig_type->index_in_settings = lig_types_and_nrs.size();
    lig_types_and_nrs.emplace_back(lig_type, n_of_lig);
}
