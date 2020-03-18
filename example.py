import matplotlib.pyplot as plt
from extrv_engine import Parameters, SimulationState, SlipBondType, CatchSlipPselBondType, CatchSlipIntegrinBondType


if __name__ == '__main__':
    p = Parameters(
        r_c=4.5,
        visc=0.01,
        temp=310,
        dens_diff=0.05,
        f_rep_0=1e3,
        tau=5
    )

    # For this parameters simulation makes sense.
    psgl = Parameters.LigandType()
    psgl_plus_esel_bond = SlipBondType(
        eq_bond_len=77,
        spring_const=100,
        binding_rate_0=0.06,
        rec_dens=3600,
        react_compl_slip=0.18,
        rup_rate_0_slip=2.6
    )
    psgl.add_bond_type(psgl_plus_esel_bond)

    # Find good parameters before you use it!
    """
    psgl_plus_psel_bond = CatchSlipPselBondType(
        eq_bond_len=77,
        spring_const=100,
        binding_rate_0=0.06,
        rec_dens=3600,
        react_compl_slip=0.18,
        rup_rate_0_slip=2.6,
        react_compl_catch = 0.18,
        rup_rate_0_catch = 2.6
    )
    """
    p.add_ligands(lig_type=psgl, n_of_lig=10000)

    # Find good parameters before you use it!
    """
    active_lfa = Parameters.LigandType()
    active_lfa_plus_icam_bond = CatchSlipIntegrinBondType(
        eq_bond_len=77,
        spring_const=100,
        binding_rate_0=0.06,
        rec_dens=3600,
        react_compl_slip=0.18,
        rup_rate_0_slip=2.6,
        react_compl_catch=0.18,
        rup_rate_0_catch=2.6
    )
    active_lfa.add_bond_type(active_lfa_plus_icam_bond)
    p.add_ligands(lig_type=active_lfa, n_of_lig=10000)
    """

    # Find good parameters before you use it!
    """
    inactive_lfa = Parameters.LigandType()
    inactive_lfa_plus_icam_bond = CatchSlipIntegrinBondType(
        eq_bond_len=77,
        spring_const=100,
        binding_rate_0=0.06,
        rec_dens=3600,
        react_compl_slip=0.18,
        rup_rate_0_slip=2.6,
        react_compl_catch=0.18,
        rup_rate_0_catch=2.6
    )
    inactive_lfa.add_bond_type(inactive_lfa_plus_icam_bond)
    p.add_ligands(lig_type=inactive_lfa, n_of_lig=10000)
    """

    s = SimulationState(h_0=0.075, p=p, seed=12345)
    # You will need more steps and smaller dt
    s.simulate(n_steps=int(1e5), dt=1e-5, shear=0.0)
    sim_hist = s.simulate_with_history(n_steps=int(5e5), dt=1e-5, shear=100)

    plt.subplot(211)
    plt.plot(sim_hist.h)
    plt.subplot(212)
    plt.plot(sim_hist.rot)
    plt.show()