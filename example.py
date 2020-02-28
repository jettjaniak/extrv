from extrv_engine import BondParameters, LigandType, ModelParameters, Settings, SimulationState


if __name__ == '__main__':
    esel_bond_p = BondParameters(
        bond_type=BondParameters.BondType.esel,
        lambda_=77,
        sigma=100,
        k_f_0=0.06,
        rec_dens=3600,
        x1s=0.18,
        k01s=2.6
    )
    psgl_lig_t = LigandType()
    psgl_lig_t.add_bond_p(esel_bond_p)

    p = ModelParameters(r_c=4.5, mu=0.01, temp=310, dens_diff=0.05, f_rep_0=1e3, tau=5)
    settings = Settings(p)
    settings.add_lig_type(psgl_lig_t, 10000)

    s = SimulationState(h_0=0.075, settings=settings, seed=12345)
    s.simulate_one_step(dt=1e-6, shear=100)