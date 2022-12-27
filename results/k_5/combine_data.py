import pandas as pd
import numpy as np
import os

ks = [5]

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)



for k in ks:
    df = pd.read_csv(os.path.join(path, f'./k_{k}/iterative_tmp.res'), sep=',|\s+', engine='python')
    merged_df = pd.DataFrame(columns=df.columns)
    Ms = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    betas = np.linspace(0.1, 1, 10)
    gammas = np.linspace(0.1, 1, 10)
    for i, beta in enumerate(betas):
        for j, gammas in enumerate
        for t in range(t_max):
            ind = i*t_max + t
            tmp = df[(df['algo'] == t) & (df['p_mask'] == p)]
            if (not tmp.empty):
                merged_df.loc[ind] = tmp.iloc[0]

                merged_df.loc[ind]['no_test'] = tmp['no_test'].sum()
                merged_df.loc[ind]['no_success'] = tmp['no_success'].sum()
                merged_df.loc[ind]['no_stop'] = tmp['no_stop'].sum()
                merged_df.loc[ind]['p_log'] = merged_df.loc[ind]['no_success']/merged_df.loc[ind]['no_test']
    
    print(len(merged_df))
    merged_df.to_csv(os.path.join(path, f'./k_{k}/combined_res.res'), index=False)