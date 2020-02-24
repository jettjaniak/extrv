from wrappers.engine import *
p = ModelParameters(4.5, 0.01, 310, 0.05, 1e3, 5)
settings = Settings(p)
psgl_lig_p = LigandType()
esel_bond_p = BondParameters("esel", 77, 100, 0.06, 3600, 0.18, 2.6)
psgl_lig_p.add_bond_p(esel_bond_p)
settings.add_lig_type(psgl_lig_p, 10000)
ss = SimulationState(0.075, settings, 12345)
sim_res = ss.simulate_with_history(1e6, 1e-6, 100, 100)