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

snps = []

for h in df.Haplotype:
    if h == '2' or h == '2-6' or h == '6':
        snps.append('A')
    elif h == '2-4' or h == '4-6':
        snps.append('A/T')
    elif h == '4':
        snps.append('T')

snps = []

for h in df.Haplotype:
    if h == '2' or h == '2-4' or h == '4':
        snps.append('T')
    elif h == '2-6' or h == '4-6':
        snps.append('T/C')
    elif h == '6':
        snps.append('C')

df['SNP 199'] = snps
df['SNP 148'] = snps

# sns.set(style='darkgrid')
# sns.boxplot(data=df, x=df['SNP 199'], y=df.Usage)
# plt.title('The distribution of usages across SNP states')
# plt.tight_layout()
# plt.show()
# plt.savefig(pth+'results/boxplot2.pdf', format='pdf')

sns.set(style='darkgrid')
sns.boxplot(data=df, x=df['SNP 148'], y=df.Usage)
plt.title('The distribution of usages across SNP states')
plt.tight_layout()
# plt.show()
plt.savefig(pth+'results/boxplot3.pdf', format='pdf')
