import matplotlib.pyplot as plt
import numpy as np
from extrv_engine import Parameters, SimulationState, SlipBondType, CatchSlipPselBondType, CatchSlipIntegrinBondType


if __name__ == '__main__':
    N_STEPS_FALLING = int(1e5)
    N_STEPS_TEST = int(5e6)

    p = Parameters(
        r_c=4.5,
        visc=0.01,
        temp=310,
        dens_diff=0.05,
        f_rep_0=1e3,
        tau=5  # TODO: use old nonspecific forces
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

    FORCES = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8]
    COEFF = 1e-5 * psgl_plus_esel_bond.binding_rate_0
    DT = 0.1 / psgl_plus_esel_bond.binding_rate_0

    seed = np.random.randint(1000)
    print("seed", seed)
    s = SimulationState(h_0=0.0242, p=p, seed=seed)
    # You will need more steps and smaller dt
    # s.simulate(n_steps=int(1e5), dt=DT, shear=0.0)
    sim_hist = s.simulate_with_history(n_steps=int(2e5), dt=DT, shear=0)#2*COEFF)
    print(len(s.bd_lig_ind))

    plt.subplot(211)
    plt.plot(sim_hist.h)
    plt.subplot(212)
    plt.plot(sim_hist.rot)
    plt.show()
