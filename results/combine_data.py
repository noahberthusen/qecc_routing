import pandas as pd
import numpy as np
import os

# type = "res_files_config"
# type = "res_files_random"
type = "res_files_dist/a_4"
ks = [5]

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)



for k in ks:
    df = pd.DataFrame()
    for file in os.listdir(os.path.join(path, f'./k_{k}/{type}')):
        if "cpp" in file:
            df = pd.concat([df, pd.read_csv(os.path.join(path, f'./k_{k}/{type}/{file}'), sep=',|\s+', engine='python')])
    merged_df = []

    Ms = np.linspace(10, 200, 20)
    betas = np.linspace(0, 1, 11)
    gammas = np.linspace(0, 1, 11)
    for m in Ms:
        for beta in betas:
            for gamma in gammas:
                tmp = df[(df['m'] == m) & (df['beta'] == round(beta,1)) & (df['gamma'] == round(gamma,1))]
                if (not tmp.empty):
                    no_test = tmp['no_test'].sum()
                    mean = sum([r['no_test'] * r['mean'] for _, r in tmp.iterrows()]) / no_test
                    variance = sum([r['no_test'] * (r['variance'] + (r['mean'] - mean)**2) for _, r in tmp.iterrows()]) / no_test

                    merged_df.append([m, k, round(beta,1), round(gamma,1), no_test, mean, variance])

    print(len(merged_df))
    merged_df = pd.DataFrame(merged_df, columns=df.columns)
    merged_df.to_csv(os.path.join(path, f'./k_{k}/{type}/combined_res.res'), index=False)