#include "helpers.h"
#include "Settings.h"
#include "AbstractBondType.h"

AbstractBondType::AbstractBondType(double lambda_, double sigma, double k_f_0, double rec_dens,
        double x1s, double k01s):

        // conversion to μm
        lambda_(lambda_ * 1e-3),
        // conversion to kg/s^2
        sigma(sigma * 1e-3),

        k_f_0(k_f_0),
        rec_dens(rec_dens),
        // conversion to μm
        x1s(x1s * 1e-4),
        k01s(k01s)
        {}

double AbstractBondType::binding_rate(double surface_dist, double temp) {
    double deviation = abs(surface_dist - lambda_);
    double rate_0 = rec_dens * k_f_0;
    // Here we assume that binding is subject to reactive compliance of slip part of bond.
    return helpers::bell_binding_rate(deviation, rate_0, sigma, x1s, temp);
}

double AbstractBondType::bond_force(double bond_length) {
    return sigma * abs(bond_length - lambda_);
}

SlipBondType::SlipBondType(double lambda_, double sigma, double k_f_0, double rec_dens, double x1s, double k01s)
        : AbstractBondType(lambda_, sigma, k_f_0, rec_dens, x1s, k01s) {}

double SlipBondType::rupture_rate(double bond_length, double temp) {
    double bond_f = bond_force(bond_length);
    // TODO: move it here from helpers
    return helpers::slip_rupture_rate(bond_f, k01s, x1s, temp);
}

CatchSlipBondType::CatchSlipBondType(double lambda_, double sigma, double k_f_0, double rec_dens, double x1s,
        double k01s, double x1c, double k01c) :

        AbstractBondType(lambda_, sigma, k_f_0, rec_dens, x1s, k01s),
        // conversion to μm
        x1c(x1c * 1e-4),
        k01c(k01c)
        {}

CatchSlipPselBondType::CatchSlipPselBondType(double lambda_, double sigma, double k_f_0, double rec_dens,
        double x1s, double k01s, double x1c, double k01c) :
        CatchSlipBondType(lambda_, sigma, k_f_0, rec_dens, x1s, k01s, x1c, k01c) {}

double CatchSlipPselBondType::rupture_rate(double bond_length, double temp) {
    double bond_f = bond_force(bond_length);
    // TODO: move it here from helpers
    return helpers::catch_slip_psel_rupture_rate(bond_f, k01s, k01c, x1s, x1c, temp);
}

CatchSlipIntegrinBondType::CatchSlipIntegrinBondType(double lambda_, double sigma, double k_f_0, double rec_dens,
        double x1s, double k01s, double x1c, double k01c) :
        CatchSlipBondType(lambda_, sigma, k_f_0, rec_dens, x1s, k01s, x1c, k01c) {}

double CatchSlipIntegrinBondType::rupture_rate(double bond_length, double temp) {
    double bond_f = bond_force(bond_length);
    // TODO: move it here from helpers
    return helpers::catch_slip_integrin_rupture_rate(bond_f, k01s, k01c, x1s, x1c, temp);
}
