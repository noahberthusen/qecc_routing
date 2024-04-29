import numpy as np
import matplotlib.pyplot as plt
import os

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

colors = [(0, 190, 255), (221, 179, 16), (0, 178, 93), (181, 29, 20), (64, 83, 211),  (251, 73, 176), ]
# colors = [(239, 230, 69), (233, 53, 161), (0, 227, 255), (225, 86, 44), (83, 126, 255), (0, 203, 133)]
# colors = [(86, 100, 26), (192, 175, 251), (230, 161, 118), (0, 103, 138), (152, 68, 100), (94, 204, 171)]
colors = [(c[0]/255, c[1]/255, c[2]/255) for c in colors]

labels = ["[[90, 8, 6]]", "[[144, 12, 12]]"]

def total_time(t, sr_time, lr_time, lr_round):
    s = 0
    for i in range(1,t+1):
        if (i % lr_round) == 0:
            s += 2*(sr_time+lr_time)
        else:
            s += 2*sr_time
    return s

def plot_savings(ax, sr_time, lr_time, i):
    t = 10000
    lr_times = np.arange(1,100,1)
    tts = [total_time(t, sr_time, lr_time, i) for i in lr_times]
    pcs = [((tts[0]-tt)/tts[0])*100 for tt in tts[1:]]

    theory_max = (2*(sr_time+lr_time) - 2*sr_time) / (2*(sr_time+lr_time)) * 100

    ax.axhline(y=theory_max, color=colors[3], linestyle='-', linewidth=1.5)
    ax.plot(lr_times[1:], pcs, color='k', linewidth=1.5, linestyle='--' if i else '-.', label=labels[i])


plt.rc('font', family='serif')
# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.linewidth'] = 1

fig, ax = plt.subplots(1, 1, figsize=(4.2,2.2), sharey=True)
plot_savings(ax, 9, 5, 0)
plot_savings(ax, 16, 15, 1)
ax.axvline(5, c='gray', alpha=0.5) #, linestyle='--')


handles,labels = ax.get_legend_handles_labels()
handles = handles[::-1]
labels = labels[::-1]
ax.legend().get_frame().set_linewidth(1)
ax.legend(handles, labels, loc='lower right', framealpha=1)

ax.set_ylabel("Circuit depth savings (%)")
ax.set_xlabel("Long-range unmasking freq. (rounds)")
plt.xscale('log')

# plt.show()
plt.savefig(os.path.join(path, f'../pot_savings.png'), dpi=1000, transparent=False, bbox_inches='tight')
