#include "Parameters.h"
#include "AbstractBondType.h"
#include "AdaptiveSimulationState.h"

#include <iostream>

int main() {
    for (int i = 26; i < 100; i++) {
        Parameters p = Parameters(
                4.5,
                0.01,
                310,
                0.05,
                5074.616349947093,
                1146.409200818891,
                1e-10,
                1e-6
        );

        Parameters::LigandType psgl = Parameters::LigandType();

        SlipBondType psgl_plus_esel_bond = SlipBondType(
                27,
                100,
                0.06,
                750,
                0.18,
                2.6
        );
        psgl.add_bond_type(&psgl_plus_esel_bond);
        p.add_ligands(&psgl, 20000);

        auto s = AdaptiveSimulationState(0.03, &p, i);
//        for (int j = 0; j < 10000; j++)
//            s.simulate_one_step();
        s.simulate(1.0, size_t(1e6));
        s.shear_rate = 0.5;
        double max_time = 4.0;
        s.simulate(max_time, size_t(3e6));
        if (s.time < max_time)
            std::cout << "only " << s.time << " seconds done!" << std::endl;
        std::cout << "seed " << i << " done" << std::endl;
//        s.simulate(5.0, size_t(5e6));
//        s.simulate(5.0, size_t(5e6));
//        s.simulate_with_history(n_steps, 0);
    }
}