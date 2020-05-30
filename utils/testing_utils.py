from collections import namedtuple
from extrv_engine import SlipBondType, Parameters
from typing import Dict, List
import numpy as np
import pandas as pd
from collections import defaultdict


SimulationStats = namedtuple("SimulationStats", ['mean_h', 'rot', 'dist', 'n_bonds', 'mean_bond_ls', 'comp_time'])
TestResultsType = Dict[float, List[SimulationStats]]


def setup_parameters(*, rec_dens=750, binding_rate_0=0.06):
    p = Parameters(
        r_cell=4.5,
        visc=0.01,
        temp=310,
        dens_diff=0.05,
        rep_0=5075,
        rep_scale=1145
    )
    psgl_plus_esel_bond = SlipBondType(
        eq_bond_len=27,
        spring_const=100,
        binding_rate_0=binding_rate_0,
        rec_dens=rec_dens,
        react_compl_slip=0.18,
        rup_rate_0_slip=2.6
    )
    psgl = Parameters.LigandType()
    psgl.add_bond_type(psgl_plus_esel_bond)
    p.add_ligands(lig_type=psgl, n_of_lig=20000)
    # Don't let them be garbage collected!
    return p, psgl, psgl_plus_esel_bond


def get_dfs_for_fields(test_results: TestResultsType, ignore=('comp_time',)):
    fields = tuple(filter(lambda f: f not in ignore, SimulationStats._fields))
    dfs_for_fields = {}
    for field in fields:
        arr = []
        columns = test_results.keys()
        for stats_list in test_results.values():
            arr.append([getattr(item_stat, field) for item_stat in stats_list])
        arr = np.array(arr).T
        dfs_for_fields[field] = pd.DataFrame(arr, columns=columns)
    return dfs_for_fields

    # arr_field_indices = tuple(enumerate(
    #     filter(lambda f_i: SimulationStats._fields[f_i] not in ignore, range(n_all_fields))
    # ))
    # n_fields = len(arr_field_indices)
    # for var, stats_list in test_results.items():
    #     n_items = len(stats_list)
    #     arr = np.empty((n_items, n_fields))
    #     for arr_i, field_i in arr_field_indices:
    #         arr[:, arr_i] = [stat[field_i] for stat in stats_list]
    #     arrays_for_fields.append(arr)
    # return np.hstack(arrays_for_vars)
