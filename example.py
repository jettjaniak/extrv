import matplotlib.pyplot as plt
import numpy as np
from extrv_engine import Parameters, SimulationState, SlipBondType, CatchSlipPselBondType, CatchSlipIntegrinBondType


if __name__ == '__main__':
    p = Parameters(
        r_cell=4.5,
        visc=0.01,
        temp=310,
        dens_diff=0.05,
        f_rep_0=1e3,
        tau=5
    )

    # For this parameters simulation makes sense.
    psgl = Parameters.LigandType()
    psgl_plus_esel_bond = SlipBondType(
        eq_bond_len=20,
        spring_const=2,
        binding_rate_0=1.2e5,
        rec_dens=1,
        react_compl_slip=0.75,
        rup_rate_0_slip=0.01
    )
    psgl.add_bond_type(psgl_plus_esel_bond)
    p.add_ligands(psgl, 10000)


    s = SimulationState(h_0=0.0225, p=p, seed=101)
    # You will need more steps and smaller dt
    s.simulate(n_steps=int(1e5), dt=0.1 / psgl_plus_esel_bond.binding_rate_0, shear=0.0)
    sim_hist = s.simulate_with_history(
        n_steps=int(5e6),
        dt=0.1 / psgl_plus_esel_bond.binding_rate_0,
        shear=3.5 * 1e-5 * psgl_plus_esel_bond.binding_rate_0,
        stop_if_no_bonds=False,
        save_every=100
    )

    plt.subplot(211)
    plt.plot(sim_hist.h)
    plt.subplot(212)
    plt.plot(sim_hist.rot)
    plt.show()
