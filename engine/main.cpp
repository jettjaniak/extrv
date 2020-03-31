#include "Parameters.h"
#include "headers/SimulationState.h"
#include "AbstractBondType.h"

#include <iostream>

int main() {
    for (int i = 0; i < 1; i++) {
        auto p = new Parameters(
                4.5,
                0.01,
                310,
                0.05,
                0.004213242356967143,
                738.979592388574
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

        std::cout << i << std::endl;
        auto s = SimulationState(0.03, p, 100 + i);
//        for (int j = 0; j < 10000; j++)
//            s.simulate_one_step();
        s.simulate_with_history(1.0, size_t(1e5), 1e-2);
        s.simulate_with_history(1.0, size_t(1e5), 1e-2);
//        s.simulate_with_history(n_steps, 0);
    }
}