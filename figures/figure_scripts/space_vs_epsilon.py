import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

colors = [(0, 190, 255), (221, 179, 16), (0, 178, 93), (181, 29, 20), (64, 83, 211),  (251, 73, 176), ]
# colors = [(239, 230, 69), (233, 53, 161), (0, 227, 255), (225, 86, 44), (83, 126, 255), (0, 203, 133)]
# colors = [(86, 100, 26), (192, 175, 251), (230, 161, 118), (0, 103, 138), (152, 68, 100), (94, 204, 171)]
colors = [(c[0]/255, c[1]/255, c[2]/255) for c in colors]

fig, ax = plt.subplots(1, 1, figsize=(4,2.75), sharey=True)

surface8x = [248,392,568,776]
surface8y = [1.4e-3,2.0e-4,9.5e-5,2.0e-5]

surface12x = [588,852,1164]
surface12y = [3.0e-4,1.4e-4,2.9e-5]

bb8x = [288,360,480,600]
bb8y = [1.6e-3,8.9e-4,1.2e-4,5.3e-5]

bb12x = [576,784]
bb12y = [1.6e-4,7.9e-5]

bb8xnoidle = [288,360,480,600]
bb8ynoidle = [3.6e-4,2.0e-4,1.3e-5,5.9e-6]

surface8xnoidle = [392,568,776,1016]
surface8ynoidle = [2.0e-4,9.5e-5,1.9e-5,7.5e-6]

hfont = {'fontname':'serif'}

plt.rcParams["font.family"] = "serif"

plt.rcParams['axes.linewidth'] = 1

s8, = ax.plot(surface8x, surface8y, marker="o", markersize=3, c=colors[0], linestyle='--',
          label="Surface $k=$8")
bb8, = ax.plot(bb8x, bb8y,  marker="o", markersize=3, c=colors[0],
         label="BB $k=$8")

s12, = ax.plot(surface12x, surface12y,  marker="o", markersize=3, c=colors[3], linestyle='--',
         label="Surface $k=$12")
bb12, = ax.plot(bb12x, bb12y,  marker="o", markersize=3, c=colors[3],
         label="BB $k=$12")

s8no, = ax.plot(surface8xnoidle, surface8ynoidle, marker="o", markersize=3, c='k', linestyle='--',
        label="Surface $k=$8 no idle")
bb8no, = ax.plot(bb8xnoidle, bb8ynoidle, marker="o", markersize=3, c='k',
        label="BB $k=$8 no idle")

ax.set_yscale('log')

leg1 = ax.legend(loc="upper right", fontsize=7, handlelength=2.7, markerscale=0.8, framealpha=0.5)
# leg2 = ax.legend([bb8no, s8no], ["BB $k=$8 no idle", "Surface $k=$8 no idle"], loc="lower right", fontsize=8, handlelength=2.7, markerscale=0.8)
ax.add_artist(leg1)
plt.xticks(fontname='serif')
plt.yticks(fontname='serif')

ax.set_ybound(3e-3,4e-6)
# ax.set_yticks([1e-3, 1e-4], **hfont)
ax.set_xlabel("Total qubits used", **hfont)
ax.set_ylabel("Logical error rate per round, $\epsilon_L$", **hfont)

plt.savefig(os.path.join(path, f'../space_vs_epsilon.png'), dpi=1000, transparent=False, bbox_inches='tight')
