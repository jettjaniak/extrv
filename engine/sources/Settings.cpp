#include "Settings.h"

#include "helpers.h"


//////////////////////
//  BondParameters  //
//////////////////////

BondParameters::BondParameters(BondType bond_type, double lambda_, double sigma, double k_f_0, double rec_dens,
                               double x1s, double k01s, double x1c, double k01c) :
        bond_type(bond_type),
        // conversion to μm
        lambda_(lambda_ * 1e-3),
        // conversion to kg/s^2
        sigma(sigma * 1e-3),

        k_f_0(k_f_0),
        rec_dens(rec_dens),
        // conversion to μm
        x1s(x1s * 1e-4),
        k01s(k01s),

        // conversion to μm
        x1c(x1c * 1e-4),
        k01c(k01c)
        {}


double BondParameters::binding_rate(double surface_dist, double temp) {
    double deviation = std::abs(surface_dist - lambda_);
    double rate_0 = rec_dens * k_f_0;
    // Here we assume that binding is subject to reactive compliance of slip part of bond.
    return helpers::bell_binding_rate(deviation, rate_0, sigma, x1s, temp);
}

double BondParameters::rupture_rate(double bond_length, double temp) {
    double bond_f = bond_force(bond_length);
    switch (bond_type) {
        case BondType::esel:
            // PSGL + E-selectin slip bond
            return helpers::esel_rupture_rate(bond_f, k01s, x1s, temp);
        case BondType::psel:
            // PSGL + P-selectin catch-slip bond
            return helpers::psel_rupture_rate(bond_f, k01s, k01c, x1s, x1c, temp);
        case BondType::integrin:
            return helpers::integrin_rupture_rate(bond_f, k01s, k01c, x1s, x1c, temp);
        default:
            abort();  // not implemented
    }
}

double BondParameters::bond_force(double bond_length) {
    return sigma * std::abs(bond_length - lambda_);
}


//////////////////
//  LigandType  //
//////////////////


vector<double> LigandType::binding_rates(double surface_dist, double temp) {
    vector<double> ret(bonds_p.size());
    BondParameters* bond_p;
    for (int i = 0; i < bonds_p.size(); i++) {
        bond_p = bonds_p[i];
        ret[i] = bond_p->binding_rate(surface_dist, temp);
    }
    return ret;
}

void LigandType::add_bond_p(BondParameters *bond_p) {
    bonds_p.push_back(bond_p);
}


///////////////////////
//  ModelParameters  //
///////////////////////

ModelParameters::ModelParameters(double r_c, double mu, double temp, double dens_diff, double f_rep_0, double tau) :
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


//////////////////////////
//  Settings  //
//////////////////////////

Settings::Settings(ModelParameters *p) : p(p) {}

void Settings::add_lig_type(LigandType *lig_type, size_t n_of_lig) {
    if (lig_type->index_in_settings > -1)
        abort();  // this ligand type was added to Settings previously, not implemented
    lig_type->index_in_settings = lig_types_and_nrs.size();
    lig_types_and_nrs.emplace_back(lig_type, n_of_lig);
}
