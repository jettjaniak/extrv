#include "../headers/Settings.h"

#include "../headers/helpers.h"


//////////////////////
//  BondParameters  //
//////////////////////

BondParameters::BondParameters(BondType bond_type_, double lambda__, double sigma_, double k_f_0_, double rec_dens_, double x1s_,
                               double k01s_, double x1c_, double k01c_) {
    bond_type = bond_type_;

    // conversion to μm
    lambda_ = lambda__ * 1e-3;
    // conversion to kg/s^2
    sigma = sigma_ * 1e-3;

    k_f_0 = k_f_0_;
    rec_dens = rec_dens_;

    // conversion to μm
    x1s = x1s_ * 1e-4;
    k01s = k01s_;

    // conversion to μm
    x1c = x1c_ * 1e-4;
    k01c = k01c_;
}

double BondParameters::binding_rate(double surface_dist, double temp) {
    double deviation = std::abs(surface_dist - lambda_);
    double rate_0 = rec_dens * k_f_0;
    // Here we assume that binding is subject to reactive compliance of slip part of bond.
    return helpers::bell_binding_rate(deviation, rate_0, sigma, x1s, temp);
}

double BondParameters::rupture_rate(double bond_length, double temp) {
    double bond_f = bond_force(bond_length);
    switch (bond_type) {
        case ESEL_BOND:
            // PSGL + E-selectin slip bond
            return helpers::esel_rupture_rate(bond_f, k01s, x1s, temp);
        case PSEL_BOND:
            // PSGL + P-selectin catch-slip bond
            return helpers::psel_rupture_rate(bond_f, k01s, k01c, x1s, x1c, temp);
        case INTEGRIN_BOND:
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

ModelParameters::ModelParameters(double r_c_, double mu_, double temp_, double dens_diff_,
        double f_rep_0_, double tau_) {
    r_c = r_c_;
    // conversion to kg/(μm s)
    mu = mu_ * 1e-7;
    temp = temp_;
    // conversion to kg/μm^3
    dens_diff = dens_diff_ * 1e-15;

    // conversion to kg μm / s^2
    f_rep_0 = f_rep_0_ * 1e-6;
    tau = tau_;
}


//////////////////////////
//  Settings  //
//////////////////////////

Settings::Settings(ModelParameters *p_) {
    p = p_;
}

void Settings::add_lig_type(LigandType *lig_type, size_t n_of_lig) {
    if (lig_type->index_in_settings > -1)
        abort();  // this ligand type was added to Settings previously, not implemented
    lig_type->index_in_settings = lig_types_and_nrs.size();
    lig_types_and_nrs.emplace_back(lig_type, n_of_lig);
}
