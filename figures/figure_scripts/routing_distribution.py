import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# from scipy.optimize import curve_fit
import lmfit
import scipy.integrate as integrate
import os

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

k = 5
type = "a_4"
df = pd.read_csv(os.path.join(path, f"../../results/k_{k}/res_files_dist/{type}/combined_res.res"))

fig, ax = plt.subplots(1, 11, figsize=(15, 1), sharex=True, sharey=False)


def fit_fun(x, c, b):
    # g(y) is the distribution function
    def g(x):
        # return 1
        a = 4
        return (a/(1-np.exp(-a))) * np.exp(-a*x)

    integrand = lambda y: g(y) * (x**y)

    return c * integrate.quad(integrand, 0, 1)[0] * x**(2*b-2)

def fit_fun_vec(x, c, b):
    if len(x) > 1:
        return [fit_fun(m, c, b) for m in x]
    else:
        return fit_fun(x, c, b)

betas = np.linspace(0, 1, 11)
gammas = [0]
for i, beta in enumerate(betas):
    for j, gamma in enumerate(gammas):
        tmp_df = df[(df['beta'] == round(beta,1)) & (df['gamma'] == round(gamma,1))]
        ax[i].scatter(tmp_df['m'], tmp_df['mean'])

        ax[i].set_xlabel('$\\beta$ = ' + f"{round(beta, 1)}")
        if (i > 7):
            ax[i].set_facecolor('xkcd:salmon')

            model = lmfit.Model(fit_fun_vec)
            params = lmfit.Parameters()
            params.add('c', value=1, min=0, max=1e10)
            params.add('b', value=beta, vary=False)

            result = model.fit(tmp_df['mean'], params, x=tmp_df['m'])
            # popt, pcov = curve_fit(fit_fun_vec, list(tmp_df['m']), list(tmp_df['mean']), bounds=([1e10,0], [beta,beta-1e-6]))
            # # ax[i].text(0.08, 0.7, round(popt[1], 2), transform=ax[i].transAxes)

            # print(beta, gamma, popt, np.sqrt(np.diag(pcov)))
            # xx = np.linspace(10, 120, 100)
            # yy = fit_fun_vec(xx, *popt)
            # result.plot_fit()
            print(result.fit_report())
            ax[i].plot(tmp_df['m'], result.best_fit, c='k')

fig.supxlabel('Grid size')
fig.supylabel('Average rounds to route')

plt.savefig(os.path.join(path, f'../routing_dist__k_{k}_{type}.png'), dpi=1000, transparent=False, bbox_inches='tight')

# plt.show()