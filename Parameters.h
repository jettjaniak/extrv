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
};

struct Parameters {
    // TODO: description and constructor with unit conversion
    double temp;
    double r_c;
    double dens_diff;
    double mu;
    
};
