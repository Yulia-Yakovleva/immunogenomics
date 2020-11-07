import pandas as pd
import itertools
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns


def get_df(vnames, jnames, c):
    df = pd.DataFrame(index=vnames, columns=jnames).fillna(0)
    for k, v in c.items():
        row_name = k.split('/')[0]
        col_name = k.split('/')[1]
        if row_name == 'IGHV1-3':
            print(k, v)
        val = v
        df.loc[row_name, col_name] = val
    print(df)
    return df


def draw_heatmap(df):
    mi = min(df.min())
    ma = max(df.max())
    sns.heatmap(df, vmin=mi, vmax=ma, cmap='rocket_r')
    plt.show()
    # plt.savefig(f"./results/{t}/hm/{prefix}_{param}.{format}", format=format)


df = pd.read_csv('analyzer_naive/cdr_details.txt', sep='\t')

v_hits = []
j_hits = []
pairs = []

for v, j in zip(df.V_hit, df.J_hit):
    v_hits.append(v.split('*')[0])
    j_hits.append(j.split('*')[0])
    pairs.append(f"{v.split('*')[0]}/{j.split('*')[0]}")

draw_heatmap(get_df(sorted(list(set(v_hits))), sorted(list(set(j_hits))), Counter(pairs)))
print(len(Counter(pairs)))