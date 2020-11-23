import time
import pandas as pd
import itertools
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO


def get_seq_pairs(names):
    c = itertools.combinations(names, 2)
    return c


def hamming_distance(seq1, seq2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))


df = pd.read_csv('/home/yulia/immunogenomics/HW2/results/analyzer_results/cdr_details.txt', sep='\t')
df['cdr3_length'] = df['CDR3_end'] - df['CDR3_start'] + 1
df['V_gene'] = [el.split('*')[0] for el in df['V_hit']]

# G = nx.Graph()

# without any optimization
# print('without optimization')
# start_time = time.time()
# for name1, seq1 in zip(df['Read_name'], df['CDR3_nucls']):
#     for name2, seq2 in zip(df['Read_name'], df['CDR3_nucls']):
#         if name1 == name2:
#             pass
#         else:
#             # 2b compute Hamming distances for pairs of computed CDR3s
#             hm = hamming_distance(seq1, seq2) / len(seq1)
#             if hm <= 0.2:
#                 # 2c connect two sequences if similarity <= 0.2
#                 G.add_edge(seq1, seq2)
# print(f"Working time is {(time.time() - start_time)}")

G = nx.Graph()

# optimization by length and drop non-unique sequences
print('with optimization')
start_time = time.time()
for length, sub_df in df.groupby('cdr3_length'):
    sub_df = sub_df.drop_duplicates(subset='CDR3_nucls').reset_index()
    G.add_nodes_from(set(sub_df['CDR3_nucls']))

    for name1, seq1 in zip(sub_df['Read_name'], sub_df['CDR3_nucls']):
        for name2, seq2 in zip(sub_df['Read_name'], sub_df['CDR3_nucls']):
            if name1 == name2:
                pass
            else:
                # 2b compute Hamming distances for pairs of computed CDR3s
                hm = hamming_distance(seq1, seq2) / length
                if hm <= 0.2:
                    # 2c connect two sequences if similarity <= 0.2
                    G.add_edge(seq1, seq2)
                    # print(name1, name2, hm)
print(f"Working time is {(time.time() - start_time)}")

# 2d report connected components
lineages = list(nx.connected_components(G))
largest = max([len(el) for el in lineages])

# 3 analyze the computed clonal lineages
print(f"Number of clonal lineages is {len(lineages)}")
print(f"Number of sequences in the larges lineage is {largest}")
print(f"The number of clonal lineages presented by at least 10 sequences is {len([el for el in lineages if len(el) <= 10])}")

# 4 usage plot of the computed V genes
v_hits = pd.DataFrame(columns=['Lineage', 'Gene'])

for num, lineage in enumerate(lineages):
    for seq in lineage:
        gene = df[df.CDR3_nucls == seq].drop_duplicates(subset='CDR3_nucls').V_gene.values[0]
        v_hits = v_hits.append({'Lineage': f"lineage_{num}", 'Gene': gene}, ignore_index=True)

counts = {}
for g in set(v_hits['Gene']):
    counts[g] = len(set(v_hits[v_hits['Gene'] == g].Lineage))

counts = pd.DataFrame.from_dict(counts, orient='index').reset_index()
counts = counts.rename(columns={"index": "Genes", 0: "Counts"}, errors="raise").sort_values(by='Counts', ascending=False)

sns.set(style='darkgrid')
sns.set_color_codes("pastel")

sns.catplot(x='Genes', y='Counts', data=counts, palette="Set3", kind='bar', height=5, aspect=15/5)
plt.xticks(rotation=90)
plt.title('Usage plot of the computed V genes')
plt.tight_layout()
plt.savefig(f"/home/yulia/immunogenomics/HW2/results/countplot.png", format='png')
plt.show()


# 5 sequences for logo
with open('/home/yulia/immunogenomics/HW2/results/seqs_for_logo.fasta', 'w') as out_f:
    for l in lineages:
        if len(l) == largest:
            reads_names = []
            for num, seq in enumerate(l):
                out_f.write(f">seq_{num}\n{seq}\n")
# 6 extract names of VDJ sequences
                for v, name in zip(df['CDR3_nucls'], df['Read_name']):
                    if seq == v:
                        reads_names.append(name)

with open('/home/yulia/immunogenomics/HW2/results/seqs_for_phylo.fasta', 'w') as out_f:
    for rec in SeqIO.parse('/home/yulia/immunogenomics/HW2/results/analyzer_results/cleaned_sequences.fasta', 'fasta'):
        if rec.id in set(reads_names):
            SeqIO.write(rec, out_f, 'fasta')