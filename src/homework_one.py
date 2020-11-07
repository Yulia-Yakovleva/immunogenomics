import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

### task1 ###

# prefs = ['naive', 'plasma', 'memory']
# prefix = 'naive'

# def get_df(vnames, jnames, c, prefix):
#     df = pd.DataFrame(index=vnames, columns=jnames).fillna(0)
#     for k, v in c.items():
#         row_name = k.split('/')[0]
#         col_name = k.split('/')[1]
#         val = v
#         df.loc[row_name, col_name] = val
#     df.to_csv(f"/home/yulia/immunogenomics/results/{prefix}_matrix.tsv", sep='\t')
#     return df
#
#
# def draw_heatmap(df, prefix):
#     mi = min(df.min())
#     ma = max(df.max())
#     sns.heatmap(df, vmin=mi, vmax=ma, cmap='rocket_r')
#     plt.title(f'VJ pairs occuring in {prefix} cells sample')
#     plt.tight_layout()
#     plt.savefig(f"/home/yulia/immunogenomics/results/{prefix}_heatmap.png", format='png')
#
#
# df = pd.read_csv(f"/home/yulia/immunogenomics/results/analyzer_{prefix}/cdr_details.txt", sep='\t')
#
# v_hits = []
# j_hits = []
# pairs = []
#
# for v, j in zip(df.V_hit, df.J_hit):
#     v_hits.append(v.split('*')[0])
#     j_hits.append(j.split('*')[0])
#     pairs.append(f"{v.split('*')[0]}/{j.split('*')[0]}")
#
# draw_heatmap(get_df(sorted(list(set(v_hits))), sorted(list(set(j_hits))), Counter(pairs), prefix=prefix), prefix=prefix)
# print(prefix, len(Counter(pairs)), max(Counter(pairs)), max(Counter(pairs).values()))

### task2 ###


# def get_idx(recs, substr_in_read):
#     for r in recs:
#         if substr_in_read in r.description:
#             return r.description.split('|')[0]
#
#
# df = pd.read_csv(f"/home/yulia/immunogenomics/results/analyzer_{prefix}/cdr_details.txt", sep='\t')
# recs = [el for el in SeqIO.parse(f"/home/yulia/immunogenomics/results/analyzer_{prefix}/v_alignment.fasta", 'fasta')]
#
# gene_to_reads = {}
#
# v_hits = []
# for v, j in zip(df.V_hit, df.J_hit):
#     v_hits.append(v.split('*')[0])
#
# for k, v in Counter(v_hits).most_common(10): # top 10
#     # print(k, v)
#     gene_to_reads[k] = []
#
# for read, gene in zip(df.Read_name, df.V_hit):
#     if gene.split('*')[0] in set(gene_to_reads.keys()):
#         gene_to_reads[gene.split('*')[0]].append(read)
#
# mut = {}
# for gene, reads in gene_to_reads.items():
#     mut[gene] = []
#     for r in reads:
#         idx = get_idx(recs, r)
#         pair = []
#         for rec in recs:
#             if rec.description.split('|')[0] == idx:
#                 pair.append(rec)
#         c = 0
#         for l1, l2 in zip(pair[0].seq, pair[1].seq):
#             if l1 != l2:
#                 c += 1
#         # print(gene, idx, round(c/len(pair[1]), 2))
#         mut[gene].append(round(c/len(pair[1]), 2))

# print(mut)
# mut_df = pd.DataFrame(columns=['V_gene', 'mutability'])
# for k, v in mut.items():
#     for el in v:
#         # print(k, el)
#         mut_df = mut_df.append({'V_gene': f"{k} n={len(v)}", 'mutability': el}, ignore_index=True)
# print(mut_df)
# mut_df.to_csv(f"/home/yulia/immunogenomics/results/{prefix}_mutability.tsv", sep='\t')
# sns.set(style='darkgrid')
# sns.boxplot(data=mut_df, x=mut_df.V_gene, y=mut_df.mutability)
# plt.xticks(rotation=45)
# plt.title(f'Mutability of top 10 V genes occuring in {prefix} cells sample')
# plt.show()
# plt.tight_layout()
# plt.savefig(f"/home/yulia/immunogenomics/results/{prefix}_mutability.png", format='png')

### task3 ###

# lens = [len(r) for r in SeqIO.parse(f"/home/yulia/immunogenomics/results/analyzer_{prefix}/cdr3s.fasta", 'fasta')]
# sns.set(style='darkgrid')
# sns.distplot(lens, bins=15)
# plt.ylabel('Number of sequences')
# plt.xlabel('Length')
# plt.title(f"CDR3 sequences length distribution of {prefix} cells")
# plt.tight_layout()
# # plt.show()
# plt.savefig(f"/home/yulia/immunogenomics/results/{prefix}_CDR3dist.png", format='png')

### task4 ###

data = pd.DataFrame(columns=['cells', 'productive'])
for prefix in ['naive', 'plasma', 'memory']:
    df = pd.read_csv(f"/home/yulia/immunogenomics/results/analyzer_{prefix}/cdr_details.txt", sep='\t')
    for val in df.Productive:
        data = data.append({'cells': prefix, 'productive': val}, ignore_index=True)

# data = data.set_index(data.columns[0])
print(data)
data = data.replace(0, 'Non-productive')
data = data.replace(1, 'Productive')

sns.set(style='darkgrid')
sns.set_color_codes("pastel")
sns.countplot(y='productive', hue='cells', data=data, palette="Set3")
plt.ylabel('Cell fraction')
plt.xlabel('Fraction')
plt.title('Fraction of productive sequences for each sample')
plt.tight_layout()
# plt.show()
plt.savefig(f"/home/yulia/immunogenomics/results/productive_fraction.png", format='png')