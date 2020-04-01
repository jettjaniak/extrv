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
                5074.616349947093,
                1146.409200818891,
                1e-10,
                1e-6
        );

        auto psgl = new Parameters::LigandType();
        double fold_change = 1.0;
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

        auto s = SimulationState(0.03, p, 110);
//        for (int j = 0; j < 10000; j++)
//            s.simulate_one_step();
        s.simulate(1.0, size_t(1e6));
        s.shear_rate = 0.5;
        s.simulate(5.0, size_t(5e6));
        s.simulate(5.0, size_t(5e6));
//        s.simulate_with_history(n_steps, 0);
    }
}