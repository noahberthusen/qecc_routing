import pandas as pd
import numpy as np
import os

codes = [
    "12_3"
]
# adv = 1.475
adv = 1.2

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

for code in codes:
    df = pd.read_csv(os.path.join(path, f'./{code}/tmp.res'), sep=',|\s+', engine='python')
    merged_df = pd.DataFrame(columns=df.columns)
    t_max = 60
    p_masks = [0.7, 0.8, 0.9, 0.95]
    # p_masks = [0.25]
    for i, p in enumerate(p_masks):
        for t in range(t_max):
            ind = i*t_max + t
            tmp = df[(df['t'] == t) & (df['p_mask'] == p) & (df['adv'] == adv)]
            if (not tmp.empty):
                merged_df.loc[ind] = tmp.iloc[0]

                merged_df.loc[ind]['no_test'] = tmp['no_test'].sum()
                merged_df.loc[ind]['no_success'] = tmp['no_success'].sum()
                merged_df.loc[ind]['p_log'] = merged_df.loc[ind]['no_success']/merged_df.loc[ind]['no_test']

    merged_df = merged_df.astype({'t':'int','k':'int','l':'int','m':'int','a1':'int','a2':'int','a3':'int','b1':'int','b2':'int','b3':'int',
                                  'emb_m':'int','emb_l':'int','aa':'int','bb':'int','adv':'float','p_phys':'float','p_mask':'float','no_test':'int','no_success':'int','p_log':'float'})
    print(len(merged_df))

    merged_df.to_csv(os.path.join(path, f'./{code}/mask_vs_synd_error_{adv}.res'), index=False)
    # merged_df.to_csv(os.path.join(path, f'./{code}/full_circuit_results.res'), index=False)
