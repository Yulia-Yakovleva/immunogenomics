import pandas as pd
import scipy.stats as stats
import itertools
import seaborn as sns
import matplotlib.pyplot as plt

pth = '/home/yulia/immunogenomics/HW3/'

df = pd.read_csv(pth+'data/IGHV1-2_usage.tsv', sep='\t')

haplotypes = sorted(list(set(df.Haplotype)))

### 1

with open(pth+'results/table1.tsv', 'w') as out_f:
    out_f.write('Haplotype\t# individuals\tMean usage\n')
    for h in haplotypes:
        subdf = df[df.Haplotype == h]
        out_f.write(f"Haplotype {h}\t{len(subdf.Haplotype)}\t{subdf.Usage.mean()}\n")

### 2

pairs = list(itertools.combinations(haplotypes, 2))

anowa_rslt = pd.DataFrame()

for p1, p2 in pairs:
    F, p = stats.f_oneway(df[df.Haplotype == p1].Usage, df[df.Haplotype == p2].Usage)
    anowa_rslt.loc[f"Haplotype {p1}", f"Haplotype {p2}"] = p
    anowa_rslt.loc[f"Haplotype {p2}", f"Haplotype {p1}"] = p

anowa_rslt.fillna('-', inplace=True)
anowa_rslt.to_csv(pth+'results/anova_rslt.tsv', sep='\t')

### 3

sns.set(style='darkgrid')
sns.boxplot(data=df, x=df.Haplotype, y=df.Usage)
plt.xticks(rotation=45)
plt.title('Usages across haplotypes')
plt.tight_layout()
# plt.show()
plt.savefig(pth+'results/boxplot.pdf', format='pdf')