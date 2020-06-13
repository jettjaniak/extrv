import os
import pickle
from collections import namedtuple
from extrv_engine import SlipBondType, Parameters
from typing import Dict, List
import numpy as np
import pandas as pd
import pingouin as pg
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy import interpolate


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


FIELDS_MAP = dict(
    dist='distance traveled in flow direction',
    rot='cell rotation',
    n_bonds='number of all formed bonds',
    mean_h='mean distance from wall'
)

FIELDS_Y_LABEL_MAP = dict(
    dist='distance in $\\mu m$',
    rot='absolute rotation in radians',
    n_bonds='number of bonds',
    mean_h='mean distance in $\\mu m$'
)

DF_COLUMNS = ("fold change", *FIELDS_MAP.keys())


def get_err_or_dt_dict(results_dir, err_or_dt, base, min_exp_bound=-float('inf'), ignore_exp=()):
    err_or_dt_dict = defaultdict(list)
    min_exp = float('inf')
    max_exp = -float('inf')
    for tr_fn in os.listdir(results_dir):
        tr_path = os.path.join(results_dir, tr_fn)
        with open(tr_path, 'rb') as tr_file:
            tr = pickle.load(tr_file)
        slice_pos = 3 if err_or_dt == 'err' else 2
        dash_split = tr_fn.split('_', maxsplit=1)[0]
        if '^' not in dash_split:
            continue
        base_str, exp_str = dash_split.split('^')
        tr_base = float(base_str[slice_pos:])
        if tr_base != base:
            continue
        tr_exp = float(exp_str)
        if tr_exp < min_exp_bound or tr_exp in ignore_exp:
            continue
        min_exp = min(min_exp, tr_exp)
        max_exp = max(max_exp, tr_exp)
        for fold_change, stats_l in tr.items():
            for stat in stats_l:
                tr_row = [fold_change] + [abs(getattr(stat, field)) for field in FIELDS_MAP.keys()]
                err_or_dt_dict[tr_exp].append(tr_row)

    err_or_dt_dict = {exp: pd.DataFrame(arr, columns=DF_COLUMNS) for exp, arr in sorted(err_or_dt_dict.items())}
    print(err_or_dt_dict.keys())
    return err_or_dt_dict, min_exp, max_exp


def get_svd(err_or_dt, base, var, min_exp_bound=-float('inf'), ignore_exp=()):
    err_or_dt_dict, min_exp, max_exp = get_err_or_dt_dict(f'results/diff_dens/diff_{err_or_dt}', err_or_dt, base,
                                                          min_exp_bound=min_exp_bound, ignore_exp=ignore_exp)

    sorted_steps = sorted(err_or_dt_dict.keys())
    steps_dfs = {s: err_or_dt_dict[s] for s in sorted_steps}
    steps_var_data = {s: sdf.groupby('fold change')[var].apply(np.array) for s, sdf in steps_dfs.items()}
    # mapping from dt / err to DF with fold change in columns
    steps_var_data = {
        s: pd.DataFrame({i: svd[i] for i in svd.index}) for s, svd in steps_var_data.items()
    }

    return steps_var_data


def eq_means_test(err_or_dt, var, min_exp_bound=-float('inf'), ignore_exp=()):
    steps_var_data = get_svd(err_or_dt, var, min_exp_bound, ignore_exp)
    svd_items = tuple(steps_var_data.items())
    print("Normality")
    for s, df in svd_items:
        print(s)
        print(pg.multivariate_normality(df))
    print("T-test")
    for (s1, df1), (s2, df2) in zip(svd_items[:-1], svd_items[1:]):
        print(s1, s2)
        print(pg.multivariate_ttest(df1, df2, paired=True))
        print()

    return steps_var_data


def mwu_test(err_or_dt, base, var, fold_change, min_exp_bound=-float('inf'), ignore_exp=()):
    steps_var_data = get_svd(err_or_dt, base, var, min_exp_bound, ignore_exp)
    svd_items = tuple(steps_var_data.items())
    print("Mann-Whitney U Test ")
    for (s1, df1), (s2, df2) in zip(svd_items[:-1], svd_items[1:]):
        print(s1, s2)
        print(pg.mwu(df1[fold_change], df2[fold_change]))
        print()


def wilcoxon_test(err_or_dt, base, min_exp_bound=-float('inf'), ignore_exp=()):
    dict_for_steps_pairs = defaultdict(lambda: defaultdict(list))
    for var in FIELDS_MAP.keys():
        steps_var_data = get_svd(err_or_dt, base, var, min_exp_bound, ignore_exp)
        svd_items = tuple(steps_var_data.items())
        for (s1, df1), (s2, df2) in zip(svd_items[:-1], svd_items[1:]):
            dict_for_this_pair = dict_for_steps_pairs[(s1, s2)]
            for fc in range(1, 6):
                w_test_res = pg.wilcoxon(df1[fc], df2[fc])
                pval = w_test_res['p-val'].values[0]
                dict_for_this_pair[var].append(pval)
    return {k: pd.DataFrame(v, index=range(1, 6)) for k, v in dict_for_steps_pairs.items()}


def kw_test(err_or_dt, var, min_exp_bound=-float('inf'), ignore_exp=()):
    err_or_dt_dict, _min_exp, _max_exp = get_err_or_dt_dict(f'results/diff_dens/diff_{err_or_dt}', err_or_dt,
                                                            min_exp_bound=min_exp_bound, ignore_exp=ignore_exp)
    step_items = tuple(err_or_dt_dict.items())
    for s, df in step_items:
        df['step'] = np.repeat(s, len(df))
    for (s1, df1), (s2, df2) in zip(step_items[:-1], step_items[1:]):
        print(s1, s2)
        data = df1.append(df2)
        print(pg.kruskal(data, dv=var, between='step'))
        print()


def compare_hists(err_or_dt, var, min_exp_bound=-float('inf'), ignore_exp=()):
    svd = get_svd(err_or_dt, var, min_exp_bound, ignore_exp)
    for fc in range(1, 6):
        plt.figure()
        plt.title(f"fold change {fc}")
        for (step, df), color in zip(svd.items(), plt.get_cmap('tab10').colors):
            hist, edges = np.histogram(df[fc])
            x = edges[:-1] + np.diff(edges)
            transparent_color = list(color[:3]) + [0.15]
            x_ls = np.linspace(np.min(x), np.max(x), 300)
            spline = interpolate.pchip(x, hist)
            hist_smooth = spline(x_ls)
            plt.plot(x_ls, hist_smooth, color=color, label=step)
            plt.fill_between(x_ls, hist_smooth, color=transparent_color)
        plt.legend()

