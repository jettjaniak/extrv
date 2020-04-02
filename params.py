from extrv_engine import SlipBondType, Parameters

REP_0 = 5074.616349947093
REP_SCALE = 1146.409200818891


def setup_parameters(*, fold_change=1.0, rep_0=REP_0, rep_scale=REP_SCALE, abs_err=1e-10, rel_err=1e-6):
    p = Parameters(
        r_cell=4.5,
        visc=0.01,
        temp=310,
        dens_diff=0.05,
        rep_0=rep_0,
        rep_scale=rep_scale,
        abs_err=abs_err,
        rel_err=rel_err
    )
    psgl_plus_esel_bond = SlipBondType(
        eq_bond_len=27,
        spring_const=100,
        binding_rate_0=0.06,
        rec_dens=750 * fold_change,
        react_compl_slip=0.18,
        rup_rate_0_slip=2.6
    )
    psgl = Parameters.LigandType()
    psgl.add_bond_type(psgl_plus_esel_bond)
    p.add_ligands(lig_type=psgl, n_of_lig=20000)
    # Don't let them be garbage collected!
    return p, psgl, psgl_plus_esel_bond