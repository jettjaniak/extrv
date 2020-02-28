#include "headers/Settings.h"
#include "headers/SimulationState.h"

#include <iostream>

int main() {
    auto p = new ModelParameters(4.5, 0.01, 310, 0.05, 1e3, 5);
    auto settings = new Settings(p);

    auto psgl_lig_t = new LigandType();
    auto psgl2_lig_t = new LigandType();
    auto esel_bond_p = new BondParameters(BondParameters::BondType::esel, 77, 100, 0.06, 3600, 0.18, 2.6);
    auto esel2_bond_p = new BondParameters(BondParameters::BondType::esel, 77, 100, 0.06, 3600, 0.18, 2.6);
    psgl_lig_t->add_bond_p(esel_bond_p);
    psgl2_lig_t->add_bond_p(esel2_bond_p);

    settings->add_lig_type(psgl_lig_t, 20000);
    settings->add_lig_type(psgl2_lig_t, 20000);

    auto s = SimulationState(0.0745478, settings, 1234567);
    s.simulate_with_history(1e5, 1e-6, 0);
}