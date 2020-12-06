import pandas as pd
import scipy.stats as stats
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
pd.set_option('display.max_columns',10)

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

anowa_rslt.fillna(1, inplace=True)
anowa_rslt.to_csv(pth+'results/anova_rslt.tsv', sep='\t')

### 3

# sns.set(style='darkgrid')
# sns.boxplot(data=df, x=df.Haplotype, y=df.Usage)
# plt.title('Usages across haplotypes')
# plt.tight_layout()
# # plt.show()
# plt.savefig(pth+'results/boxplot.pdf', format='pdf')

# видим, что четвертый гаплотип отличается от всех
# print(anowa_rslt < 0.05)

### 5

df = df.set_index('SubjectID')
# print(df)

IGHV4_T199 = df[(df.Haplotype == '4') | (df.Haplotype == '2-4') | (df.Haplotype == '4-6')]
IGHV4_T199 = IGHV4_T199.replace('4', 'T199').replace('2-4', 'T199').replace('4-6', 'T199')

IGHV6_C148 = df[(df.Haplotype == '6') | (df.Haplotype == '2-6') | (df.Haplotype == '4-6')]
IGHV6_C148 = IGHV6_C148.replace('6', 'C148').replace('2-6', 'C148').replace('4-6', 'C148')

snps = pd.concat([IGHV4_T199, IGHV6_C148])

sns.set(style='darkgrid')
sns.boxplot(data=snps, x=snps.Haplotype, y=snps.Usage)
plt.title('The distribution of usages across SNP states')
plt.xlabel('SNPs')
plt.tight_layout()
# plt.show()
plt.savefig(pth+'results/boxplot2.pdf', format='pdf')
