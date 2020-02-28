#include "headers/Settings.h"
#include "headers/SimulationState.h"

#include <iostream>

int main() {
    auto p = new ModelParameters(4.5, 0.01, 310, 0.05, 1e3, 5);
    auto settings = new Settings(p);

    auto psgl_lig_t = new LigandType();
    auto psgl2_lig_t = new LigandType();
    auto esel_bond_p = new BondParameters(BondParameters::BondType::esel, 77, 100, 0.06, 3600, 0.18, 2.6);
    auto esel2_bond_p = new BondParameters(BondParameters::BondType::esel, 77, 100, 0.06, 3600, 0.18, 2.6);
    psgl_lig_t->add_bond_p(esel_bond_p);
    psgl2_lig_t->add_bond_p(esel2_bond_p);

    settings->add_lig_type(psgl_lig_t, 20000);
    settings->add_lig_type(psgl2_lig_t, 20000);

    auto s = SimulationState(0.0745478, settings, 1234567);
    size_t n_steps = 1e5;

    auto hist = History(&s);
    for (size_t i = 0; i < n_steps; ++i) {
        s.simulate_one_step(1e-6, 0);
        if (i % (n_steps/100) == 0) {
            hist.update();
            std::cout << "h: " << s.h << std::endl;
            std::cout << "rot: " << s.rot * (180 / PI) << std::endl;
            std::cout << std::endl;
        }
    }
    hist.finish();
    std::cout << "bonds trajectories:" << std::endl;
    for (auto &final_traj : hist.final_trajs_vec) {
        std::cout << final_traj.start_i << ", " << final_traj.n_of_pos << std::endl;
    }
}