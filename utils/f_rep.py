import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

G = 9806650
PI = 3.141592653589793238463
g = (4 / 3) * PI * 4.5**3 * 0.05 * 1e-15 * G

print(g)

def function_factory(stat_h, f):
    def function(tau):
        return np.log(1 + f * tau / g) - tau * stat_h
    return function


def rep_fun(f, tau, h):
    return f * tau * np.exp(- tau * h) / (1 - np.exp(-tau * h))


def plot_f_rep(stat_h, h_min=0.02, h_max=0.04):
    ls = np.linspace(h_min, h_max, 100)
    f_min = g * stat_h * (1 + 5e-4)
    plt.figure()
    plt.yscale('log')
    for f in np.linspace(f_min, f_min * 1e12, 10):
        fun = function_factory(stat_h, f)
        tau_min = 0.01
        tau_max = 1e7

        # print(f"fun({tau_min}) = {fun(tau_min)}, fun({tau_max}) = {fun(tau_max)}")
        sol = optimize.root_scalar(fun, bracket=[tau_min, tau_max], method='brentq')
        tau = sol.root
        print(f, tau)
        plt.plot(ls, [rep_fun(f, tau, h) for h in ls], label=f"f = {int(f)}, $\\tau = {tau}$")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    plot_f_rep(0.0271)