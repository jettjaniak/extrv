#pragma once

#include "types.h"
#include "AbstractBondType.h"


struct Parameters {
    /**
     * Type of ligand. Defined by types of bonds it can form.
     */
    struct LigandType {
        /**
         * it's easier to keep pointer to parameters here, so we don't have to keep it in thousands of Ligands
         * or pass it many times in simulation
         */
        Parameters* p = nullptr;
        /// types of bonds that a ligand of this type can form
        vector<AbstractBondType*> bonds_types;
        /// maximal distance from surface s.t. rate of binding is >= MIN_RATE
        double max_surf_dist;

        /**
         * Add type of bond that a ligand of this type can form.
         */
        void add_bond_type(AbstractBondType* bond_type);

        /**
         * Compute rates of binding for each bond type that a ligand of this type can form.
         *
         * @param surface_dist ligand's distance from surface in μm
         * @return vector of binding rates for each bond type, in the same order they were added
         */
        void compute_binding_rates(double surface_dist, vector<double> &binding_rates);

        /**
         * Update maximal surface distance for which any binding rate is higher than MIN_RATE.
         */
        void update_max_surf_dist();
    };

    /// cell (sphere) radius in μm
    double r_cell;
    /// viscosity in kg/(μm s)
    double visc;
    /// temperature in K
    double temp;
    /// density difference between cell and fluid (cell is more dense) in kg/μm^3
    double dens_diff;

    /// second pair element is number of ligands of particular type on whole sphere
    vector<pair<LigandType*, size_t>> lig_types_and_nrs;

    // TODO: documentation
    double max_surf_dist = 0.0;

    /**
     * @param r_cell cell (sphere) radius in μm
     * @param visc fluid viscosity in g / (cm s)
     * @param temp temperature in K
     * @param dens_diff density difference between cell and fluid (cell is more dense) in g/cm^3
     */
    Parameters(double r_cell, double visc, double temp, double dens_diff);

    /**
     * Add n_of_lig ligands of lig_type type.
     */
    void add_ligands(LigandType *lig_type, size_t n_of_lig);
};