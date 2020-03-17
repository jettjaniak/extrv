#include "AbstractBondType.h"

#include "helpers.h"
#include "Parameters.h"


////////////////////////
//  AbstractBondType  //
////////////////////////

AbstractBondType::AbstractBondType(
        double eq_bond_len, double spring_const, double binding_rate_0,
        double rec_dens, double react_compl_slip, double rup_rate_0_slip):

        // conversion to μm
        eq_bond_len(eq_bond_len * 1e-3),
        // conversion to kg/s^2
        spring_const(spring_const * 1e-3),

        binding_rate_0(binding_rate_0),
        rec_dens(rec_dens),
        // conversion to μm
        react_compl_slip(react_compl_slip * 1e-4),
        rup_rate_0_slip(rup_rate_0_slip) {}


double AbstractBondType::binding_rate(double surface_dist, double temp) {
    double deviation = std::abs(surface_dist - eq_bond_len);
    double rate_0 = rec_dens * binding_rate_0;
    // Here we assume that binding is subject to reactive compliance of slip part of bond.
    return helpers::bell_binding_rate(deviation, rate_0, spring_const, react_compl_slip, temp);
}

double AbstractBondType::bond_force(double bond_length) {
    return spring_const * std::abs(bond_length - eq_bond_len);
}


////////////////////
//  SlipBondType  //
////////////////////

SlipBondType::SlipBondType(
        double eq_bond_len, double spring_const, double binding_rate_0,
        double rec_dens, double react_compl_slip, double rup_rate_0_slip) :

        // use base class constructor
        AbstractBondType(
                eq_bond_len, spring_const, binding_rate_0,
                rec_dens, react_compl_slip, rup_rate_0_slip) {}


double SlipBondType::rupture_rate(double bond_length, double temp) {
    double bond_f = bond_force(bond_length);
    return rup_rate_0_slip * exp((react_compl_slip * bond_f) / (K_B * temp));
}


/////////////////////////////////
//  AbstractCatchSlipBondType  //
/////////////////////////////////

AbstractCatchSlipBondType::AbstractCatchSlipBondType(
        double eq_bond_len, double spring_const, double binding_rate_0,
        double rec_dens, double react_compl_slip,  double rup_rate_0_slip,
        double react_compl_catch, double rup_rate_0_catch) :

        // use base class constructor
        AbstractBondType(
                eq_bond_len, spring_const, binding_rate_0,
                rec_dens, react_compl_slip, rup_rate_0_slip),

        // then initialize two additional parameters
        react_compl_catch(react_compl_catch * 1e-4),  // conversion to μm
        rup_rate_0_catch(rup_rate_0_catch) {}


/////////////////////////////
//  CatchSlipPselBondType  //
/////////////////////////////

CatchSlipPselBondType::CatchSlipPselBondType(
        double eq_bond_len, double spring_const, double binding_rate_0,
        double rec_dens, double react_compl_slip, double rup_rate_0_slip,
        double react_compl_catch, double rup_rate_0_catch) :

        // use base class constructor
        AbstractCatchSlipBondType(
                eq_bond_len, spring_const, binding_rate_0,
                rec_dens, react_compl_slip, rup_rate_0_slip,
                react_compl_catch, rup_rate_0_catch) {}


double CatchSlipPselBondType::rupture_rate(double bond_length, double temp) {
    // TODO: any gain from using static?
    double half_f = bond_force(bond_length) / 2;
    double slip_part = (2.0 / 3.0) * rup_rate_0_slip * exp(react_compl_slip * half_f / (K_B * temp));
    double catch_part = rup_rate_0_catch * exp(react_compl_catch * half_f / (K_B * temp));
    return slip_part + catch_part;
}


/////////////////////////////////
//  CatchSlipIntegrinBondType  //
/////////////////////////////////

CatchSlipIntegrinBondType::CatchSlipIntegrinBondType(
        double eq_bond_len, double spring_const, double binding_rate_0,
        double rec_dens, double react_compl_slip, double rup_rate_0_slip,
        double react_compl_catch, double rup_rate_0_catch) :

        AbstractCatchSlipBondType(
                eq_bond_len, spring_const, binding_rate_0,
                rec_dens, react_compl_slip, rup_rate_0_slip,
                react_compl_catch, rup_rate_0_catch) {}


double CatchSlipIntegrinBondType::rupture_rate(double bond_length, double temp) {
    double bond_f = bond_force(bond_length);
    return 1 / (
            (1 / rup_rate_0_slip) * exp(-react_compl_slip * bond_f / (K_B * temp)) +
            (1 / rup_rate_0_catch) * exp(-react_compl_catch * bond_f / (K_B * temp))
    );
}
