#include "headers/Settings.h"
#include "headers/SimulationState.h"
#include "AbstractBondType.h"

#include <iostream>

int main() {
    auto p = new Settings::ModelParameters(4.5, 0.01, 310, 0.05, 1e3, 5);
    auto settings = new Settings(p);

    auto psgl_lig_t = new Settings::LigandType();
    auto esel_bond_t = new SlipBondType(77, 100, 0.06, 3600, 0.18, 2.6);
    psgl_lig_t->add_bond_type(esel_bond_t);

    settings->add_lig_type(psgl_lig_t, 10000);

    auto s = SimulationState(0.0745478, settings, 1234567);
    s.simulate_with_history(1e4, 1e-5, 0);
}