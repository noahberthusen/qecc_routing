import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

def plot_surface_data_10x(ax, d, k):
    ts = [10,15,20,25,30,35,40,45,50,55,60,65,70,75,80]

    data = {
        4 : [0.0016812,0.0025489,0.0033865,0.0042501,0.005075872203551685,0.005890609331292553,0.006775329147990811,0.007623126039106398,0.008449990198795465,0.009316634479898236,0.010100351027599614,0.010910404010550826,0.011730542378952338,0.012582637703544886,0.013414575945573883],
        5: [0.0002328,0.0003646,0.0004877,0.0006244110640902478,0.000727458675529053,0.0008700268201796874,0.0010078452954527433,0.0011384220747020423,0.0012852489506495737,0.0014111965349456826,0.0015131452842902181,0.0016341289794200478,0.0017465836669805232,0.0018875359930499269,0.0020079426162348275],
        6: [0.000115,0.0001742,0.0002347,0.0002891,0.000357,0.0004147,0.0004773,0.0005367770074793773,0.0005985242143951054,0.0006724841370966704,0.0007164854556718286,0.0007670300696188646,0.0008337658125202725,0.0008905651288221894,0.0009531508929132619],
        7: [2.06e-05,3.22e-05,4.41e-05,6.07e-05,7.44e-05,9.24e-05,9.6e-05,0.0001048,0.0001281,0.0001325,0.000148,0.0001623,0.0001691,0.0001841,0.0001947]
    }
    labels  = {
        4: 'h',
        5: 's',
        6: 'D',
        7: '^'
    }
    ax.plot(ts, [1-(1-a)**k for a in data[d]], labels[d], c='k', markersize=4, label=f"Surface code: [[{d**2*k},{k},{d}]] ({(2*d**2-1)*k} qubits)")
    popt, pcov = curve_fit(fun, ts, [1-(1-a)**k for a in data[d]], maxfev=1000, p0=(0.001))
    print(d, k, popt)
    xx = np.linspace(2, 80, 1000)
    yy = fun(xx, *popt)
    ax.plot(xx, yy, c='k', linewidth=1)

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


x10 = True
colors = [(0, 190, 255), (221, 179, 16), (0, 178, 93), (181, 29, 20), (64, 83, 211),  (251, 73, 176), ]
# colors = [(239, 230, 69), (233, 53, 161), (0, 227, 255), (225, 86, 44), (83, 126, 255), (0, 203, 133)]
# colors = [(86, 100, 26), (192, 175, 251), (230, 161, 118), (0, 103, 138), (152, 68, 100), (94, 204, 171)]
colors = [(c[0]/255, c[1]/255, c[2]/255) for c in colors]

plt.rc('font', family='serif')
# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.linewidth'] = 1

fig, ax = plt.subplots(1, 3, figsize=(12,3), sharey=True)

codes = [[12,3],[9,5],[12,5],[15,5]]#,[21,5]]
ds = [6,6,8,8]
# codes = [[14,7]]

for i, code in enumerate(codes):
    if x10:
        df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/10xfull_circuit_results_5.res'))
    else:
        df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/full_circuit_results_5.res'))
    df['p_error'] = 1 - df['p_log']
    df['p_std_dev'] = np.sqrt(df['p_error'] * df['p_log'] / df['no_test'])
    # df['p_std_dev'].replace(to_replace=0, value=1e-2, inplace=True)
    guesses = []
    params = []

    def fun(x, a):
        return 1 - (1 - a)**x

    # colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    tmp_df = df[(df['p_std_dev'] > 0) & (df['p_phys'] == 0.0002) & (df['k'] == 8)]

    ax[0].errorbar(tmp_df['t'], tmp_df['p_error'], tmp_df['p_std_dev'], fmt='o', c=colors[i], markersize=5, label=f"Quasi-cyclic: [[{code[0]*code[1]*2},8,{ds[i]}]] ({code[0]*code[1]*8} qubits)")
    popt, pcov = curve_fit(fun, tmp_df['t'], tmp_df['p_error'], maxfev=1000, p0=(0.001), sigma=tmp_df['p_std_dev'])
    print(code, popt)
    xx = np.linspace(2, 80, 1000)
    yy = fun(xx, *popt)
    ax[0].plot(xx, yy, c=colors[i], linewidth=1)

if x10:
    plot_surface_data_10x(ax[0], 4, 8)
    plot_surface_data_10x(ax[0], 5, 8)
    plot_surface_data_10x(ax[0], 6, 8)
    plot_surface_data_10x(ax[0], 7, 8)
else:
    plot_surface_data(ax[0], 4, 8)
    plot_surface_data(ax[0], 5, 8)
    plot_surface_data(ax[0], 6, 8)
    plot_surface_data(ax[0], 7, 8)


ax[0].set_yscale('log')
ax[0].set_ylabel('Logical error rate, $p_\log$')
ax[0].set_xlabel('Rounds, $t$')
ax[0].legend(loc='upper center', bbox_to_anchor=(0.5,1.73), frameon=False, fontsize=9)


codes = [[12,6],[14,7]]#,[18,6]]
arr = [3,6]
ds = [12,8]

for i, code in enumerate(codes):
    if x10:
        df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/10xfull_circuit_results_5.res'))
    else:
        df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/full_circuit_results_5.res'))
    df['p_error'] = 1 - df['p_log']
    df['p_std_dev'] = np.sqrt(df['p_error'] * df['p_log'] / df['no_test'])
    # df['p_std_dev'].replace(to_replace=0, value=1e-2, inplace=True)
    guesses = []
    params = []

    def fun(x, a):
        return 1 - (1 - a)**x

    tmp_df = df[(df['p_std_dev'] > 0) & (df['p_phys'] == 0.0002) & (df['k'] == 12) & (df['a1'] == arr[i])]

    ax[1].errorbar(tmp_df['t'], tmp_df['p_error'], tmp_df['p_std_dev'], c=colors[i], markersize=5, fmt='o', label=f"Quasi-cyclic: [[{code[0]*code[1]*2},12,{ds[i]}]] ({code[0]*code[1]*8} qubits)")
    popt, pcov = curve_fit(fun, tmp_df['t'], tmp_df['p_error'], maxfev=1000, p0=(0.001), sigma=tmp_df['p_std_dev'])
    print(code, popt)
    xx = np.linspace(2, 80, 1000)
    yy = fun(xx, *popt)
    ax[1].plot(xx, yy, c=colors[i], linewidth=1)

if x10:
    plot_surface_data_10x(ax[1], 5, 12)
    plot_surface_data_10x(ax[1], 6, 12)
    plot_surface_data_10x(ax[1], 7, 12)
else:
    plot_surface_data(ax[1], 5, 12)
    plot_surface_data(ax[1], 6, 12)
    # plot_surface_data(ax[1], 7, 12)



ax[1].set_yscale('log')
ax[1].set_xlabel('Rounds, $t$')
ax[1].legend(loc='upper center', bbox_to_anchor=(0.5,1.5), frameon=False, fontsize=9)

# handles, labels = plt.gca().get_legend_handles_labels()
# order = [3,2,1,0]
# plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])


# codes = [[15,5]]

# for code in codes:
#     if x10:
#         df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/10xfull_circuit_results_5.res'))
#     else:
#         df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/full_circuit_results_5.res'))
#     df['p_error'] = 1 - df['p_log']
#     df['p_std_dev'] = np.sqrt(df['p_error'] * df['p_log'] / df['no_test'])
#     # df['p_std_dev'].replace(to_replace=0, value=1e-2, inplace=True)
#     guesses = []
#     params = []

#     def fun(x, a):
#         return 1 - (1 - a)**x

#     colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
#     tmp_df = df[(df['p_std_dev'] > 0) & (df['p_phys'] == 0.0002) & (df['k'] == 16)]

#         # tmp_df_fit = df[(df['p_mask'] == j) & (df['algo'] >= 200)]
#         # tmp_df_before = df[(df['p_mask'] == j) & (df['algo'] < 200) & (df['algo'] > 10)]

#     ax[2].errorbar(tmp_df['t'], tmp_df['p_error'], tmp_df['p_std_dev'], fmt='o', label=f"Quasi-cyclic: [[{code[0]*code[1]*2},16,d]] ({code[0]*code[1]*8} qubits)")
# # ax.plot(ts, p_error, '-o', c='k', label="Surface code: [[712,8,6]]")

# if x10:
#     plot_surface_data_10x(ax[2], 5, 16)
#     plot_surface_data_10x(ax[2], 6, 16)
#     plot_surface_data_10x(ax[2], 7, 16)
# else:
#     plot_surface_data(ax[2], 5, 16)
#     plot_surface_data(ax[2], 6, 16)
#     plot_surface_data(ax[2], 7, 16)


# codes = [[12,5],[12,5],[9,5],[9,5]]
codes = [[9,5],[9,5],[14,7],[14,7]]
arr = [8,8,6,6]
ds = [8,8,8,8]
ks = [8,8,12,12]

lr_round = [1,5,1,5]

for i, code in enumerate(codes):
    if x10:
        df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/10xfull_circuit_results_{lr_round[i]}.res'))
    else:
        df = pd.read_csv(os.path.join(path, f'../../src_py/quasi-cyclic/results/{code[0]}_{code[1]}/full_circuit_results_{lr_round[i]}.res'))

    df['p_error'] = 1 - df['p_log']
    df['p_std_dev'] = np.sqrt(df['p_error'] * df['p_log'] / df['no_test'])
    # df['p_std_dev'].replace(to_replace=0, value=1e-2, inplace=True)
    guesses = []
    params = []

    def fun(x, a):
        return 1 - (1 - a)**x

    tmp_df = df[(df['p_std_dev'] > 0) & (df['p_phys'] == 0.0002) & (df['a1'] == arr[i])]


    ax[2].errorbar(tmp_df['t'], tmp_df['p_error'], tmp_df['p_std_dev'], fmt='o', c=colors[i], markersize=5, label=f"Quasi-cyclic: [[{code[0]*code[1]*2},{ks[i]},{ds[i]}]] (Every {lr_round[i]} round(s))")
    popt, pcov = curve_fit(fun, tmp_df['t'], tmp_df['p_error'], maxfev=1000, p0=(0.001), sigma=tmp_df['p_std_dev'])
    print(code, popt)
    xx = np.linspace(2, 80, 1000)
    yy = fun(xx, *popt)
    ax[2].plot(xx, yy, c=colors[i], linewidth=1)

ax[2].set_yscale('log')
ax[2].set_xlabel('Rounds, $t$')
ax[2].legend(loc='upper center', bbox_to_anchor=(0.5,1.42), frameon=False, fontsize=9)

# https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot

for i, label in enumerate(('(a)', '(b)', '(c)')):
    # ax = fig.add_subplot(2,2,i+1)
    ax[i].text(-0.05, 1.08, label, transform=ax[i].transAxes,
      fontsize=14, va='top', ha='right')

plt.savefig(os.path.join(path, f'../sim_results.png'), dpi=1000, transparent=False, bbox_inches='tight')
