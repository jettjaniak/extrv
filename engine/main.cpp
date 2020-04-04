#include "Parameters.h"
#include "AbstractBondType.h"
#include "AdaptiveSimulationState.h"
#include "RKSimulationState.h"
#include "EulerSimulationState.h"

#include <iostream>

int main() {
    for (int i = 100; i < 200; i++) {
        Parameters p = Parameters(
                4.5,
                0.01,
                310,
                0.05,
                5074.616349947093,
                1146.409200818891
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

        auto s = EulerSimulationState(0.03, &p, i, 1e-5);
        size_t max_steps_falling = 2e5;
        size_t max_steps_rolling = 6e5;

//        auto s = RKSimulationState(0.03, &p, i, 1e-5);
//        size_t max_steps_falling = 2e5;
//        size_t max_steps_rolling = 6e5;

//        auto s = AdaptiveSimulationState(0.03, &p, i);
//        size_t max_steps_falling = 1e6;
//        size_t max_steps_rolling = 3e6;

        s.simulate(1.0, max_steps_falling);
        s.shear_rate = 0.5;
        double max_time = 4.0;
        s.simulate(max_time, max_steps_rolling);
        if (s.time < max_time)
            std::cout << "only " << s.time << " seconds done!" << std::endl;
        std::cout << "seed " << i << " done" << std::endl;
    }
}