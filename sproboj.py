from utils.e_sel_mean_v import iteration_wrapper
for i in range(10000):
    print(i)
    fc, mean_v = iteration_wrapper(0.714285714285713, max_time=0.1, seed=100 + i)