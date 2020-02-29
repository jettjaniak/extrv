#pragma once

#include "types.h"

struct AbstractBondType {
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
    AbstractBondType(double lambda_, double sigma, double k_f_0, double rec_dens, double x1s, double k01s);


    // TODO: documentation
    double binding_rate(double surface_dist, double temp);

    // pure virtual
    // TODO: documentation
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

    SlipBondType(double lambda_, double sigma, double k_f_0, double rec_dens, double x1s, double k01s);
};

struct CatchSlipBondType : AbstractBondType {
    // reactive compliance for catch bond in μm
    double x1c;
    // multiplicative constant in catch part of rupture rate in 1/s
    double k01c;

    CatchSlipBondType(double lambda_, double sigma, double k_f_0, double rec_dens, double x1s, double k01s,
            double x1c, double k01c);
};

struct CatchSlipPselBondType : CatchSlipBondType {
    double rupture_rate(double bond_length, double temp) override;

    CatchSlipPselBondType(double lambda_, double sigma, double k_f_0, double rec_dens, double x1s, double k01s,
            double x1c, double k01c);
};

struct CatchSlipIntegrinBondType : CatchSlipBondType {
    double rupture_rate(double bond_length, double temp) override;

    CatchSlipIntegrinBondType(double lambda_, double sigma, double k_f_0, double rec_dens, double x1s, double k01s,
            double x1c, double k01c);
};