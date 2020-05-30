#include "Parameters.h"
#include "AbstractBondType.h"
#include "EulerSS.h"
#include "RKSS.h"
#include "AdapSS.h"

#include <iostream>

int main() {
    for (int i = 220; i < 226; i++) {
        Parameters p = Parameters(
                4.5,
                0.01,
                310,
                0.05,
                5075,
                1145
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
        unsigned int seed = i;
        double dt = 1e-5;
        auto s = EulerProbSS(0.03, &p, i, 1e-5);
        double falling_time = 1.0;
        size_t max_steps_falling = falling_time / dt + 10;
        double rolling_time = 5.0;
        size_t max_steps_rolling = rolling_time / dt + 10;
        double max_time = falling_time + rolling_time;

        s.simulate(falling_time, max_steps_falling);
        std::cout << s.bd_lig_ind.size() << " bonds" << std::endl;
        s.shear_rate = 100;

        s.simulate(max_time, max_steps_rolling);
        if (s.time < max_time)
            std::cout << "only " << s.time << " seconds done!" << std::endl;
        std::cout << "seed " << seed << " done" << std::endl;
    }
}