#include "Parameters.h"
#include "headers/SimulationState.h"
#include "AbstractBondType.h"

#include <iostream>

int main() {
    for (int i = 20; i < 100; i++) {
        auto p = new Parameters(
                4.5,
                0.01,
                310,
                0.05,
                // not using it
                1e3,
                5
        );

        auto psgl = new Parameters::LigandType();
        double fold_change = 0.714285714285713;
        double REC_DENS_0 = 750;
        double BINDING_RATE_0 = 0.06;
        auto psgl_plus_esel_bond = new SlipBondType(
                27,
                100,
                BINDING_RATE_0,
                REC_DENS_0 * fold_change,
                0.18,
                2.6
        );
        psgl->add_bond_type(psgl_plus_esel_bond);
        p->add_ligands(psgl, 20000);

        double k_on_0 = fold_change * REC_DENS_0 * BINDING_RATE_0;
        double dt = std::min(0.1 / k_on_0, 1e-6);
        size_t n_steps = 0.1 / dt;
        size_t n_steps_falling = 0.1 / dt;

        std::cout << i << std::endl;
        auto s = SimulationState(0.03, p, 100 + i);
        s.rot = - PI / 2;
        s.simulate_with_history(n_steps_falling, dt, 0);
        s.simulate_with_history(n_steps, dt, 1);
    }
}