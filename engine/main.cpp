#include "SimulationSettings.h"
#include "SimulationState.h"

#include <iostream>

int main() {
    auto p = new Parameters(4.5, 0.01, 310, 0.05, 1e3, 5);
    auto settings = new SimulationSettings(p);

    auto psgl_lig_p = new LigandParameters(psgl);
    auto esel_bond_p = new BondParameters(77, 100, 0.06, 3600, 0.18, 2.6);
    psgl_lig_p->add_bond_p(esel_bond_p);

    settings->add_lig_type(psgl_lig_p, 10000);

    auto ss = SimulationState(0.0745478, settings, 1234567);

    size_t n_steps = 1e6;
    for (size_t i = 0; i < n_steps; ++i) {
        ss.simulate_one_step(1e-6, 0);
        if (i % (n_steps/10) == 0) {
            std::cout << "h: " << ss.h << std::endl;
            std::cout << "alpha_0: " << ss.alpha_0 * (180 / PI) << std::endl;
            std::cout << "bonds: " << ss.bd_lig_ind.size() << std::endl;
            std::cout << std::endl;
        }
    }
}