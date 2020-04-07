import numpy as np
import pickle
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams.update({
    'errorbar.capsize': 4,
    'font.size': 15
})


def standard_error_of_mean(values):
    return np.sqrt(np.var(values) / len(values))


def plot(all_tests_results):
    for name, test_results in all_tests_results.items():
        x = []
        y = []
        err = []
        for shear, v_list in test_results.items():
            x.append(shear)
            y.append(np.mean(v_list))
            err.append(standard_error_of_mean(v_list))

        x = np.array(x)
        y = np.array(y)
        err = np.array(err)

        sort_ind = np.argsort(x)
        x = x[sort_ind]
        y = y[sort_ind]
        err = err[sort_ind]

        plt.errorbar(x, y, err, label=name)
    plt.xlabel("shear rate $(1/s)$")
    plt.ylabel("mean velocity $(\\mu m/s)$")
    plt.tight_layout()
    plt.legend()