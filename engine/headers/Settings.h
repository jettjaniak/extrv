#pragma once

#include "types.h"


struct BondParameters {
    enum class BondType {psel, esel, integrin};  // TODO: change esel to slip, psel to catch-slip-psel, etc.

    // TODO: documentation
    BondType bond_type;
    // optimal bond length in μm
    double lambda_;
    // spring constant in kg/s^2
    double sigma;

    // multiplicative constant in binding rate in μm^2/s, it doesn't include receptor density
    double k_f_0;
    // density of bond receptors on the surface in 1/μm^2
    double rec_dens;

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
     * @param lambda_ equilibrium bond length in nm
     * @param sigma spring constant in dyn / cm
     * @param k_f_0 multiplicative constant in binding rate in μm^2/s
     * @param rec_dens density of bond receptors on the surface in 1/μm^2
     * @param x1s reactive compliance for slip bond in Å
     * @param k01s multiplicative constant in slip part of binding and rupture rate in 1/s
     * @param x1c reactive compliance for catch bond in Å
     * @param k01c multiplicative constant in catch part of rupture rate in 1/s
     */
    BondParameters(BondType bond_type, double lambda_, double sigma, double k_f_0, double rec_dens,
                   double x1s, double k01s, double x1c = 0.0, double k01c = 0.0);

    // TODO: documentation
    double binding_rate(double surface_dist, double temp);

    // TODO: documentation
    double rupture_rate(double bond_length, double temp);

    /**
     * Compute force exerted on bond.
     *
     * @param bond_length bond length in μm
     * @return force exerted on bond in kg μm / s^2
     */
    double bond_force(double bond_length);
};


struct LigandType {
    int index_in_settings = -1;
    vector<BondParameters*> bonds_p;

    LigandType() = default;

    // TODO: documentation
    void add_bond_p(BondParameters* bond_p);

    // TODO: documentation
    vector<double> binding_rates(double surface_dist, double temp);
};


struct ModelParameters {
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

    // TODO: documentation
    ModelParameters(double r_c, double mu, double temp, double dens_diff, double f_rep_0, double tau);
};


struct Settings {
    ModelParameters* p = nullptr;
    // second pair element is number of ligands of particular type on whole sphere
    vector<pair<LigandType*, size_t>> lig_types_and_nrs;

    explicit Settings(ModelParameters* p);

    // TODO: documentation
    void add_lig_type(LigandType* lig_type, size_t n_of_lig);
};