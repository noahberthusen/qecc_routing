import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import matplotlib.ticker as mticker

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)


# arr = np.array([(4, 0.49577464788732395, 0.004295749734674281, 0.98935),
#  (6, 0.5016077170418006, 0.0037279699324495913, 0.98445),
#  (8, 0.5072815533980582, 0.003726771492750664, 0.9794),
#  (10, 0.49300155520995337, 0.003890052141068276, 0.97428),
#  (12, 0.5042428198433421, 0.003992324832879425, 0.96936),
#  (14, 0.502093217973765, 0.0040656730659530994, 0.96417),
#  (16, 0.517198404785643, 0.00450056257032129, 0.95988),
#  (18, 0.4924645390070922, 0.004346095844504022, 0.95488),
#  (20, 0.4986528098537336, 0.004177038943504493, 0.94804),
#  (22, 0.5050993022007515, 0.004671065871561577, 0.94411),
#  (24, 0.49563642351391407, 0.004578023358565695, 0.93927),
#  (26, 0.4938013442867812, 0.004629976957290606, 0.93305),
#  (28, 0.4945722970039079, 0.005650385107045794, 0.93091),
#  (30, 0.4895625581704561, 0.005276873668616659, 0.92479),
#  (32, 0.4973828673560577, 0.00639057363264509, 0.92167),
#  (34, 0.4972806810120596, 0.005625832950995171, 0.91542),
#  (36, 0.5002244165170556, 0.00625768487616371, 0.91088),
#  (38, 0.5004819535182606, 0.007301765880237804, 0.90663),
#  (40, 0.5038402457757296, 0.006793372859755084, 0.90235),
#  (42, 0.49624654286843145, 0.007343450976901509, 0.89876),
#  (44, 0.49911668991166896, 0.00719368031822511, 0.89245),
#  (46, 0.5071866703286643, 0.007858369725069323, 0.89077),
#  (48, 0.5073355003074761, 0.00845210286965255, 0.88617),
#  (50, 0.5022462562396006, 0.008092748351898158, 0.8798)])

arr = np.array([(4, 0.5130784708249497, 0.0019597399073394773, 0.99503),
 (6, 0.50920245398773, 0.0019357765791198267, 0.99185),
 (8, 0.5161616161616162, 0.0018684981315018685, 0.9901),
 (10, 0.4988066825775656, 0.0019241870309794113, 0.98743),
 (12, 0.49796195652173914, 0.0019080870412471582, 0.98528),
 (14, 0.5030826140567201, 0.0016467096302018743, 0.98378),
 (16, 0.49229979466119095, 0.0018561579570024068, 0.98052),
 (18, 0.4668181818181818, 0.001768916155419223, 0.978),
 (20, 0.5119396732299958, 0.0019976847346152666, 0.97613),
 (22, 0.5030257186081695, 0.002146760343481655, 0.97356),
 (24, 0.4928893513701006, 0.0020799654025556802, 0.97117),
 (26, 0.48778135048231513, 0.0020022706161626587, 0.9689),
 (28, 0.4894376673609045, 0.0021626879417212513, 0.96639),
 (30, 0.4966298193583176, 0.002139348433394606, 0.96291),
 (32, 0.4876974876974877, 0.0024339758058644257, 0.96139),
 (34, 0.48296151017406225, 0.002293554070537213, 0.95921),
 (36, 0.49331489165514064, 0.002487926240304405, 0.95662),
 (38, 0.4945078612965755, 0.0021603028618769468, 0.95357),
 (40, 0.5015125324114088, 0.0025479176278152917, 0.95372),
 (42, 0.4905698641249239, 0.0027243370604508304, 0.95069),
 (44, 0.5027037466203167, 0.0026997954061293793, 0.94822),
 (46, 0.4981167608286252, 0.0025768296546625834, 0.9469),
 (48, 0.5007233273056058, 0.0024346353339684554, 0.9447),
 (50, 0.5013750429700928, 0.002972967233653989, 0.94182)])

plt.rc('font', family='serif')
# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.linewidth'] = 1

fig, ax1 = plt.subplots(figsize=(5,2.2))



ax2 = ax1.twinx()
ax2.ticklabel_format(style='sci', axis='y', scilimits=(-3,-3))
# f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
# g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
# plt.gca().yaxis.set_major_formatter(mticker.FuncFormatter(g))


ax1.scatter(arr[:,0], arr[:,3], c='k', marker='.', label="Success prob")
ax2.scatter(arr[:,0], arr[:,2], c='r', marker='.', label="Bell fidelity")

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# ax1.legend(loc='lower right')
ax2.legend(lines+lines2, labels+labels2, loc='upper center')

ax1.set_xlabel("Bell chain length")
ax1.set_ylabel("Success probability")
ax2.set_ylabel("Effective CNOT error rate")

plt.savefig(os.path.join(path, f'../bell_fidelity.png'), dpi=1000, transparent=False, bbox_inches='tight')
