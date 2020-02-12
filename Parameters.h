#pragma once

#include "types.h"

struct BondParameters {
    // optimal bond length in μm
    double lambda_;
    // spring constant in kg/s^2
    double sigma;

    // reactive compliance for slip bond in μm
    double x1s;
    // multiplicative constant in slip part of binding or rupture rate in 1/s
    double k01s;

    // reactive compliance for catch bond in μm
    double x1c;
    // multiplicative constant in catch part of rupture rate in 1/s
    double k01c;

    /**
     * Bond parameters constructor with unit conversion.
     *
     * @param lambda__ equilibrium bond length in nm
     * @param sigma_ spring constant in dyn / cm
     * @param x1s_ reactive compliance for slip bond in Å
     * @param k01s_ multiplicative constant in slip part of binding and rupture rate in 1/s
     * @param x1c_ reactive compliance for catch bond in Å
     * @param k01c_ multiplicative constant in catch part of rupture rate in 1/s
     */
    BondParameters(double lambda__, double sigma_, double x1s_, double k01s_, double x1c_ = 0.0, double k01c_ = 0.0) {
        // conversion to μm
        lambda_ = lambda__ * 1e-3;
        // conversion to kg/s^2
        sigma = sigma_ * 1e-3;

        // conversion to μm
        x1s = x1s_ * 1e-4;
        k01s = k01s_;

        // conversion to μm
        x1c = x1c_ * 1e-4;
        k01c = k01c_;
    }

};

struct LigandParameters {
    LigandType lig_type;
    vector<BondParameters*> bonds_p;

    explicit LigandParameters(LigandType lig_type_) {
        lig_type = lig_type_;
    }
    
    void add_bond_p(BondParameters* bond_p) {
        bonds_p.push_back(bond_p);
    }
};

struct Parameters {
    // cell radius in μm
    double r_c;
    // viscosity in kg/(μm s)
    double mu;
    // temperature in K
    double temp;
    // density difference in kg/μm^3
    double dens_diff;

    // repulsive force coefficient in kg μm / s^2
    double f_rep_0;
    // reciprocal length scale of repulsive force in Å
    double tau;

    Parameters(double r_c_, double mu_, double temp_, double dens_diff_, double f_rep_0_, double tau_) {
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
};
