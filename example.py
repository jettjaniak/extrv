from wrappers.engine import *
import matplotlib.pyplot as plt


if __name__ == '__main__':
    p = ModelParameters(4.5, 0.01, 310, 0.05, 1e3, 5)
    settings = Settings(p)
    psgl_lig_p = LigandType()
    esel_bond_p = BondParameters("esel", 77, 100, 0.06, 3600, 0.18, 2.6)
    psgl_lig_p.add_bond_p(esel_bond_p)
    settings.add_lig_type(psgl_lig_p, 20000)

    # add seed parameter (positive integer) for reproducibility
    ss = SimulationState(h_0=0.075, settings=settings)
    sim_res = ss.simulate_with_history(n_steps=1e6, dt=1e-6, shear=0, save_every=100)

    plt.subplot(211)
    plt.plot(sim_res.rot)
    plt.title("rotation")
    plt.subplot(212)
    plt.plot(sim_res.h)
    plt.title("height")

    # TODO: trajectory plot with individual bonds

    plt.tight_layout()
    plt.show()
