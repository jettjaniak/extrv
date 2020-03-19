#include "Parameters.h"
#include "headers/SimulationState.h"
#include "AbstractBondType.h"

#include <iostream>

int main() {
    auto p = new Parameters(
            4.5,
            0.01,
            310,
            0.05,
            1e3,
            5
        );

    auto psgl = new Parameters::LigandType();
    auto psgl_plus_esel_bond = new SlipBondType(
            20,
            2,
            1.2e5,
            1,
            0.75,
            0.01
        );
    psgl->add_bond_type(psgl_plus_esel_bond);

    p->add_ligands(psgl, 1);

    auto s = SimulationState(0.0242, p, 788);
    std::cout << s.ligands.size() << std::endl;
//    s.simulate(size_t(1e5), 1e-5, 0, false);
    s.simulate_with_history(size_t(2e4), 0.1 / psgl_plus_esel_bond->binding_rate_0, 0, false, 100);
}