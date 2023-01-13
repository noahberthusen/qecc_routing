import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os

k = 5
# type = "res_files_config"
type = "res_files_random"

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

df = pd.read_csv(os.path.join(path, f"../../results/k_{k}/{type}/combined_res.res"))

fig, ax = plt.subplots(10, 11, figsize=(15, 15), sharex=True, sharey=False)

betas = np.linspace(0.1, 1, 10)
gammas = np.linspace(0, 1, 11)
for i, beta in enumerate(betas):
    for j, gamma in enumerate(gammas):
        tmp_df = df[(df['beta'] == round(beta,1)) & (df['gamma'] == round(gamma,1))]
        ax[9-i][j].scatter(tmp_df['m'], tmp_df['mean'])

        if (i == 0):
            ax[9-i][j].set_xlabel('$\gamma$ = ' + f"{round(gamma, 1)}")
        if (j == 0):
            ax[9-i][j].set_ylabel('$\\beta$ = ' + f"{round(beta, 1)}")
        if (2*beta + gamma >= 2):
            ax[9-i][j].set_facecolor('xkcd:salmon')

            def fun(x, c, e):
                return c*x**e
            
            popt, pcov = curve_fit(fun, tmp_df['m'], tmp_df['mean'])
            ax[9-i][j].text(0.08, 0.7, round(popt[1], 2), transform=ax[9-i][j].transAxes)
            # ax[10-i][j].text(0.4, 0.2, round(2*beta+gamma-2,2), transform=ax[10-i][j].transAxes)

            print(beta, gamma, popt, np.sqrt(np.diag(pcov))[1])
            xx = np.linspace(10, 120, 100)
            yy = fun(xx, *popt)
            ax[9-i][j].plot(xx, yy, c='k')
        # ax[i][j].set_yscale('log')
fig.supxlabel('Grid size')
fig.supylabel('Average rounds to route')

plt.savefig(os.path.join(path, f'../routing_dist__k_{k}_{type}.png'), dpi=1000, transparent=False, bbox_inches='tight')
