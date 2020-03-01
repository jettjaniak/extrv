#pragma once

#include "types.h"
#include "AbstractBondType.h"


struct Settings {

    struct ModelParameters {
        /// cell (sphere) radius in μm
        double r_cell;
        /// viscosity in kg/(μm s)
        double visc;
        /// temperature in K
        double temp;
        /// density difference between cell and fluid (cell is more dense) in kg/μm^3
        double dens_diff;

        /// repulsive force coefficient in kg μm / s^2
        double f_rep_0;

        /**
         * Reciprocal length scale of repulsive force in Å.
         * We keep original unit because the formula for repulsive force is constructed for angstroms.
         */
        double tau;

        /**
         *
         * @param r_cell cell (sphere) radius in μm
         * @param visc fluid viscosity in g / (cm s)
         * @param temp temperature in K
         * @param dens_diff density difference between cell and fluid (cell is more dense) in g/cm^3
         * @param f_rep_0 repulsive force coefficient in N
         * @param tau reciprocal length scale of repulsive force in Å
         */
        ModelParameters(double r_cell, double visc, double temp, double dens_diff, double f_rep_0, double tau);
    };

    struct LigandType {
        // TODO: documentation
        int index_in_settings = -1;
        /// types of bonds that a ligan
        vector<AbstractBondType*> bonds_types;

        LigandType() = default;

        /**
         * Add type of bond that ligands of this type can form.
         */
        void add_bond_type(AbstractBondType* bond_type);

        ///
        vector<double> binding_rates(double surface_dist, double temp);
    };

    /// model parameters
    ModelParameters* p;
    /// second pair element is number of ligands of particular type on whole sphere
    vector<pair<LigandType*, size_t>> lig_types_and_nrs;

    explicit Settings(ModelParameters* p);

    // TODO: documentation
    void add_ligands(LigandType *lig_type, size_t n_of_lig);
};