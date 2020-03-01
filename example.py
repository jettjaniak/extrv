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

    p.add_ligands(lig_type=psgl, n_of_lig=10000)

    s = SimulationState(h_0=0.075, p=p, seed=12345)
    sim_hist = s.simulate_with_history(n_steps=int(1e5), dt=1e-6, shear=0.0)
