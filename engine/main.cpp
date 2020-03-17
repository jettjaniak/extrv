#include "Parameters.h"
#include "headers/SimulationState.h"
#include "AbstractBondType.h"

#include <iostream>

int main() {
    auto p = new Parameters(4.5, 0.01, 310, 0.05, 1e3, 5);

    auto psgl_lig_t = new Parameters::LigandType();
    auto esel_bond_t = new SlipBondType(77, 100, 0.06, 3600, 0.18, 2.6);
    psgl_lig_t->add_bond_type(esel_bond_t);

    p->add_ligands(psgl_lig_t, 10000);

    auto s = SimulationState(0.0745478, p, 1234567);
    s.simulate_with_history(size_t(1e5), 1e-5, 0);
}