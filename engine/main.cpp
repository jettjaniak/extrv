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
            5074.616349947093,
            1146.409200818891
    );

    auto psgl = new Parameters::LigandType();
    double REC_DENS_0 = 750;
    double BINDING_RATE_0 = 0.06;
    auto psgl_plus_esel_bond = new SlipBondType(
            27,
            100,
            BINDING_RATE_0,
            REC_DENS_0,
            0.18,
            2.6
    );
    psgl->add_bond_type(psgl_plus_esel_bond);
    p->add_ligands(psgl, 20000);

    double dt = 1e-6;
    size_t n_steps = 5e6;
    size_t n_steps_falling = 1e6;

    auto s = SimulationState(0.03, p, 104);
    s.simulate(n_steps_falling, dt, 0);
    std::cout << s.bd_lig_ind.size() << std::endl;
}