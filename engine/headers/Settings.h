#pragma once

#include "types.h"
#include "AbstractBondType.h"


struct Settings {

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

    struct LigandType {
        int index_in_settings = -1;
        vector<AbstractBondType*> bonds_types;

        LigandType() = default;

        // TODO: documentation
        void add_bond_type(AbstractBondType* bond_type);

        // TODO: documentation
        vector<double> binding_rates(double surface_dist, double temp);
    };

    ModelParameters* p = nullptr;
    // second pair element is number of ligands of particular type on whole sphere
    vector<pair<LigandType*, size_t>> lig_types_and_nrs;

    explicit Settings(ModelParameters* p);

    // TODO: documentation
    void add_lig_type(LigandType* lig_type, size_t n_of_lig);
};