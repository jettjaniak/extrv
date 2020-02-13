#include "SimulationSettings.h"
#include "SimulationState.h"

#include <iostream>

int main() {
    auto p = new Parameters(4.5, 0.01, 310, 0.05, 1e3, 5);
    auto settings = new SimulationSettings(p);

    auto psgl_lig_p = new LigandParameters(psgl);
    auto esel_bond_p = new BondParameters(77, 100, 3600, 0.06, 0.18, 2.6);
    psgl_lig_p->add_bond_p(esel_bond_p);

    settings->add_lig_type(psgl_lig_p, 10000);

    auto ss = SimulationState(0.9, settings, 12345);

    for (int i = 0; i < 10000; ++i) {
        ss.simulate_one_step(1e-5, 0.0);
        if (i % 1000 == 0) std::cout << ss.h << std::endl;
    }
}