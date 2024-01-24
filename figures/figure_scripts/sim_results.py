import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

def plot_surface_data(ax, d, k):
    ts = [10,20,30,40,50,60,70,80,90,100]
    data = {
        4: [0.0018403432365179412,0.003548884836582971,0.005237650883477318,0.007009329456665087,0.00840560329087456,0.0105202064820794,0.01219697646294924,0.013865183674463015,0.015314132520496883,0.01715576234983606],
        5: [0.0002723,0.0005522729795094148,0.0008386979167081266,0.0011151446076454652,0.0014134666568070217,0.0016878039413872328,0.0019517033589881271,0.002296813870443364,0.0025723677714511403,0.0028566649941461782],
        6: [0.0001253,0.0002542,0.0003836,0.0005061233422585084,0.0006370462446168446,0.0007557720022841242,0.000895625954188077,0.0010187308860719292,0.0011239277522523494,0.0012853476839334776],
        7: [2.57e-05,5.98e-05,8.27e-05,0.000113,0.0001406,0.0001667,0.0001895,0.0002301,0.0002509,0.0002893]
    }
    labels  = {
        4: '-o',
        5: '-.o',
        6: ':o',
        7: '--o'
    }
    print([1-(1-a)**k for a in data[d]][:8])
    ax.plot(ts[:8], [1-(1-a)**k for a in data[d]][:8], labels[d], c='k', label=f"Surface code: [[{d**2*k},{k},{d}]] ({(2*d**2-1)*k} qubits)")



plt.rc('font', family='serif')
# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.linewidth'] = 1

fig, ax = plt.subplots(1, 3, figsize=(12,4), sharey=True)

codes = [[12,3],[9,5],[12,5],[15,5]]#,[21,5]]
# codes = [[14,7]]

for code in codes:
    df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/full_circuit_results_5.res'))
    df['p_error'] = 1 - df['p_log']
    df['p_std_dev'] = np.sqrt(df['p_error'] * df['p_log'] / df['no_test'])
    # df['p_std_dev'].replace(to_replace=0, value=1e-2, inplace=True)
    guesses = []
    params = []

    def fun(x, a):
        return 1 - (1 - a)**x

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    tmp_df = df[(df['p_std_dev'] > 0) & (df['p_phys'] == 0.0001) & (df['k'] == 8)]

        # tmp_df_fit = df[(df['p_mask'] == j) & (df['algo'] >= 200)]
        # tmp_df_before = df[(df['p_mask'] == j) & (df['algo'] < 200) & (df['algo'] > 10)]

    ax[0].errorbar(tmp_df['t'], tmp_df['p_error'], tmp_df['p_std_dev'], fmt='o', label=f"Quasi-cyclic: [[{code[0]*code[1]*2},8,d]] ({code[0]*code[1]*8} qubits)")
# ax.plot(ts, p_error, '-o', c='k', label="Surface code: [[712,8,6]]")
plot_surface_data(ax[0], 4, 8)
plot_surface_data(ax[0], 5, 8)
plot_surface_data(ax[0], 6, 8)
plot_surface_data(ax[0], 7, 8)


# popt, pcov = curve_fit(fun, tmp_df['t'], tmp_df['p_error'], maxfev=1000, p0=(0.001),
#     sigma=tmp_df['p_std_dev'])
# xx = np.linspace(1, 100, 1000)
# yy = fun(xx, *popt)
# ax.plot(xx, yy, c='k')


# ax.plot(np.linspace(0, 0.05, 100), np.linspace(1e-3, 50*1e-3, 100), c='k')
# ax[1].plot(np.linspace(1e-3,1e-2,100), np.linspace(1e-3, 1e-2, 100), c='k')

# ax.set_title('ISD with $p_0 = 0.001$')
# ax[1].set_title('SSF with $k=1$')
# ax.legend(loc='lower right')
ax[0].set_yscale('log')
ax[0].set_ylabel('Logical error rate, $p_\log$')
ax[0].set_xlabel('Rounds, $t$')
ax[0].legend(loc='upper center', bbox_to_anchor=(0.5,1.6), frameon=False)



# codes = [[12,3],[9,5],[12,5],[15,5]]#,[21,5]]
codes = [[14,7]]

for code in codes:
    df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/full_circuit_results_5.res'))
    df['p_error'] = 1 - df['p_log']
    df['p_std_dev'] = np.sqrt(df['p_error'] * df['p_log'] / df['no_test'])
    # df['p_std_dev'].replace(to_replace=0, value=1e-2, inplace=True)
    guesses = []
    params = []

    def fun(x, a):
        return 1 - (1 - a)**x

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    tmp_df = df[(df['p_std_dev'] > 0) & (df['p_phys'] == 0.0001) & (df['k'] == 12)]

        # tmp_df_fit = df[(df['p_mask'] == j) & (df['algo'] >= 200)]
        # tmp_df_before = df[(df['p_mask'] == j) & (df['algo'] < 200) & (df['algo'] > 10)]

    ax[1].errorbar(tmp_df['t'], tmp_df['p_error'], tmp_df['p_std_dev'], fmt='o', label=f"Quasi-cyclic: [[{code[0]*code[1]*2},12,d]] ({code[0]*code[1]*8} qubits)")
# ax.plot(ts, p_error, '-o', c='k', label="Surface code: [[712,8,6]]")
# plot_surface_data(ax, 4, 12)
plot_surface_data(ax[1], 5, 12)
plot_surface_data(ax[1], 6, 12)
plot_surface_data(ax[1], 7, 12)


# popt, pcov = curve_fit(fun, tmp_df['t'], tmp_df['p_error'], maxfev=1000, p0=(0.001),
#     sigma=tmp_df['p_std_dev'])
# xx = np.linspace(1, 100, 1000)
# yy = fun(xx, *popt)
# ax.plot(xx, yy, c='k')


# ax.plot(np.linspace(0, 0.05, 100), np.linspace(1e-3, 50*1e-3, 100), c='k')
# ax[1].plot(np.linspace(1e-3,1e-2,100), np.linspace(1e-3, 1e-2, 100), c='k')

# ax.set_title('ISD with $p_0 = 0.001$')
# ax[1].set_title('SSF with $k=1$')
# ax.legend(loc='lower right')
ax[1].set_yscale('log')
# ax.set_ylabel('Logical error rate, $p_\log$')
ax[1].set_xlabel('Rounds, $t$')
ax[1].legend(loc='upper center', bbox_to_anchor=(0.5,1.55), frameon=False)
# ax[1].set_xlabel('$p$')

# handles, labels = plt.gca().get_legend_handles_labels()
# order = [3,2,1,0]
# plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])


codes = [[15,5]]

for code in codes:
    df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/full_circuit_results_5.res'))
    df['p_error'] = 1 - df['p_log']
    df['p_std_dev'] = np.sqrt(df['p_error'] * df['p_log'] / df['no_test'])
    # df['p_std_dev'].replace(to_replace=0, value=1e-2, inplace=True)
    guesses = []
    params = []

    def fun(x, a):
        return 1 - (1 - a)**x

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    tmp_df = df[(df['p_std_dev'] > 0) & (df['p_phys'] == 0.0001) & (df['k'] == 16)]

        # tmp_df_fit = df[(df['p_mask'] == j) & (df['algo'] >= 200)]
        # tmp_df_before = df[(df['p_mask'] == j) & (df['algo'] < 200) & (df['algo'] > 10)]

    ax[2].errorbar(tmp_df['t'], tmp_df['p_error'], tmp_df['p_std_dev'], fmt='o', label=f"Quasi-cyclic: [[{code[0]*code[1]*2},16,d]] ({code[0]*code[1]*8} qubits)")
# ax.plot(ts, p_error, '-o', c='k', label="Surface code: [[712,8,6]]")
# plot_surface_data(ax, 4, 12)
plot_surface_data(ax[2], 5, 16)
plot_surface_data(ax[2], 6, 16)
plot_surface_data(ax[2], 7, 16)


# popt, pcov = curve_fit(fun, tmp_df['t'], tmp_df['p_error'], maxfev=1000, p0=(0.001),
#     sigma=tmp_df['p_std_dev'])
# xx = np.linspace(1, 100, 1000)
# yy = fun(xx, *popt)
# ax.plot(xx, yy, c='k')


# ax.plot(np.linspace(0, 0.05, 100), np.linspace(1e-3, 50*1e-3, 100), c='k')
# ax[1].plot(np.linspace(1e-3,1e-2,100), np.linspace(1e-3, 1e-2, 100), c='k')

# ax.set_title('ISD with $p_0 = 0.001$')
# ax[1].set_title('SSF with $k=1$')
# ax.legend(loc='lower right')
ax[2].set_yscale('log')
# ax.set_ylabel('Logical error rate, $p_\log$')
ax[2].set_xlabel('Rounds, $t$')
ax[2].legend(loc='upper center', bbox_to_anchor=(0.5,1.55), frameon=False)

# https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot


plt.savefig(os.path.join(path, f'../sim_results.png'), dpi=1000, transparent=False, bbox_inches='tight')
