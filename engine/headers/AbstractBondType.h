#pragma once

#include "types.h"

struct AbstractBondType {
    // TODO: add name

    /// equilibrium bond length in μm
    double eq_bond_len;
    /// spring constant in kg/s^2
    double spring_const;

    /// multiplicative constant in binding rate in μm^2/s, it doesn't include receptor density
    double binding_rate_0;
    /// density of bond receptors on the surface in 1/μm^2
    double rec_dens;

    /// reactive compliance for binding rate and slip part of rupture rate in μm
    double react_compl_slip;
    /// multiplicative constant in slip part of rupture rate in 1/s
    double rup_rate_0_slip;

    /**
     * Abstract bond type constructor with unit conversion.
     *
     * @param eq_bond_len equilibrium bond length in nm
     * @param spring_const spring constant in dyn / cm
     * @param binding_rate_0 multiplicative constant in binding rate in μm^2/s, it doesn't include receptor density
     * @param rec_dens density of bond receptors on the surface in 1/μm^2
     * @param react_compl_slip reactive compliance for binding rate and slip part of rupture rate in Å
     * @param rup_rate_0_slip multiplicative constant in slip part of rupture rate in 1/s
     */
    AbstractBondType(double eq_bond_len, double spring_const, double binding_rate_0, double rec_dens,
            double react_compl_slip, double rup_rate_0_slip);


    /**
     * Compute rate of binding of this type.
     *
     * @param surface_dist distance from surface in μm
     * @param temp temperature in K
     * @return binding rate in 1/s
     */
    double binding_rate(double surface_dist, double temp);

    /**
     * Compute rate of rupture of bond of this type.
     *
     * @param bond_length bond length in μm
     * @param temp temperature in K
     * @return rupture rate in 1/s
     */
    virtual double rupture_rate(double bond_length, double temp) = 0;

    /**
     * Compute force exerted on bond.
     *
     * @param bond_length bond length in μm
     * @return force exerted on bond in kg μm / s^2
     */
    double bond_force(double bond_length);
};

struct SlipBondType : AbstractBondType {
    double rupture_rate(double bond_length, double temp) override;

    /**
     * Slip bond type constructor with unit conversion.
     * Inherited from AbstractBondType.
     */
    SlipBondType(double eq_bond_len, double spring_const, double binding_rate_0, double rec_dens,
            double react_compl_slip, double rup_rate_0_slip);
};

struct AbstractCatchSlipBondType : AbstractBondType {
    /// reactive compliance for catch bond in μm
    double react_compl_catch;
    /// multiplicative constant in catch part of rupture rate in 1/s
    double rup_rate_0_catch;

    /**
     * Abstract catch-slip bond type constructor with unit conversion.
     * All but last two parameters are passed to AbstractBondType constructor.
     *
     * @param react_compl_catch reactive compliance for catch part of rupture rate in Å
     * @param rup_rate_0_catch multiplicative constant in catch part of rupture rate in 1/s
     */
    AbstractCatchSlipBondType(double eq_bond_len, double spring_const, double binding_rate_0, double rec_dens,
            double react_compl_slip, double rup_rate_0_slip, double react_compl_catch, double rup_rate_0_catch);
};

struct CatchSlipPselBondType : AbstractCatchSlipBondType {
    double rupture_rate(double bond_length, double temp) override;

    /**
     * Catch-slip, P-selectin like bond type constructor with unit conversion.
     * Inherited from AbstractCatchSlipBondType.
     */
    CatchSlipPselBondType(double eq_bond_len, double spring_const, double binding_rate_0, double rec_dens,
            double react_compl_slip, double rup_rate_0_slip, double react_compl_catch, double rup_rate_0_catch);
};

struct CatchSlipIntegrinBondType : AbstractCatchSlipBondType {
    double rupture_rate(double bond_length, double temp) override;

    /**
     * Catch-slip, integrin like bond type constructor with unit conversion.
     * Inherited from AbstractCatchSlipBondType.
     */
    CatchSlipIntegrinBondType(double eq_bond_len, double spring_const, double binding_rate_0, double rec_dens,
            double react_compl_slip, double rup_rate_0_slip, double react_compl_catch, double rup_rate_0_catch);
};