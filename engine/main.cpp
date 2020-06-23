#include "Parameters.h"
#include "AbstractBondType.h"
#include "EulerSS.h"
#include "RKSS.h"
#include "AdapSS.h"

#include <iostream>

int main() {
    array<uint32_t, 80> seeds {1731712957, 3190037783, 3759517182, 4177028560, 1830405957,
                           3762189420, 1546919209, 752713666, 2277760481, 2225307367, 771879583, 2686016408, 3091179815,
                           54641801, 822526847, 3485949959, 1808968316, 2740946722, 2128909877, 1735928740, 709544941,
                           256163407, 2950738251, 67318337, 3705137511, 23055023, 3076940865, 2962143592, 3516257955,
                           3987088566, 997047223, 3296440446, 3213720315, 501563541, 3609636494, 3947717566, 1205502233,
                           3141042261, 2198693213, 1654518378, 8356930, 682264254, 1174718341, 556715953, 3943414791,
                           3560480401, 1483454036, 1038552105, 325198726, 712250723, 783014756, 2310304458, 1895812771,
                           2484051420, 622190755, 1641609963, 3466489807, 1441286071, 2419453959, 1533929822, 2464008272,
                           4191568628, 1751682214, 2288410354, 3596124070, 1820690884, 2414023370, 2264240200, 3133142593,
                           3102573477, 2513463839, 1460214160, 629244395, 2998528036, 1651710757, 156828413, 2688925137,
                           1635904256, 1031466595, 2716957902};
    for (auto seed : seeds) {
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
        double dt = pow(2, -23);
        auto s = EulerProbSS(0.03, &p, seed, dt);
        double falling_time = 1.0;
        size_t max_steps_falling = falling_time / dt + 10;
        double rolling_time = 10.0;
        size_t max_steps_rolling = rolling_time / dt + 10;
        double max_time = falling_time + rolling_time;

        s.simulate(falling_time, max_steps_falling);
        std::cout << s.bd_lig_ind.size() << " bonds" << std::endl;

        s.shear_rate = 10;
        s.simulate(max_time, max_steps_rolling);
        if (s.time < max_time)
            std::cout << "only " << s.time << " seconds done!" << std::endl;
        std::cout << "seed " << seed << " done" << std::endl;
    }
}