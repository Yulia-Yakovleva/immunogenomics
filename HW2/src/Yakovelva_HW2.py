import time
import pandas as pd
import itertools
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt


def get_seq_pairs(names):
    c = itertools.combinations(names, 2)
    return c


def hamming_distance(seq1, seq2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))


df = pd.read_csv('/home/yulia/immunogenomics/HW2/results/analyzer_results/cdr_details.txt', sep='\t')
df['cdr3_length'] = df['CDR3_end'] - df['CDR3_start'] + 1

G = nx.Graph()

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
                    # 2c connect two seqiences if similarity <= 0.2
                    G.add_edge(seq1, seq2)
                    # print(name1, name2, hm)
print(f"Working time is {(time.time() - start_time)}")

# 2d report connected components
lineages = list(nx.connected_components(G))

print(f"Number of clonal lineages is {len(lineages)}")
print(f"Number of sequences in the larges lineage is {max([len(el) for el in lineages])}")
print(f"The number of clonal lineages presented by at least 10 sequences is {len([el for el in lineages if len(el) <= 10])}")

        # print(i, seq1, length)
        # for j in range(i+1, length):
        #     print(sub_df['CDR3_nucls'][j])

        # print(i)
        # print(sub_df['CDR3_nucls'][i+1])
        # print(seq)
    # break
    # for i in range(length):
    #     sub_df['CDR3_nucls'][i]
    # for sequences in pairs:
    #     print(p)
    # for seq1, seq2 in pairs:
    #     print(seq1, seq2)dfgdg
    # print(pairs)
    # break

    # for

# nx.draw(G)
# plt.show()




# (2, 3, {'weight': 3.1415})


# def split_on_kmers(rec):
#     k = round(len(rec) * 0.2) + 1
#     return set(rec.seq[i: j] for i in range(len(rec.seq)) for j in range(i + 1, len(rec.seq) + 1) if len(rec.seq[i:j]) == k)


# combinations = get_seq_pairs(recs)

# start_time = time.time()
# print('Calculate with optimization...')
# for pair in combinations:
#     if len(pair[0]) == len(pair[1]):
#         kmers1 = split_on_kmers(pair[0])
#         kmers2 = split_on_kmers(pair[1])
#         if len(kmers1.intersection(kmers2)) > 0:
#             # hd = hamming_distance(pair)
#             print(pair)
# print('Done')
# print(f"Working time is {(time.time() - start_time)}")

# start_time = time.time()
# print('Calculate without optimization...')
# for pair in combinations:
#     if len(pair[0]) == len(pair[1]):
#         print(pair)
#         # hd = hamming_distance(pair)
# print('Done')
# print(f"Working time is {(time.time() - start_time)}")
#         print()
# hd = hamming_distance(pair)
# similarity = hd / len(pair[0])
# if similarity <= 0.2:
#     print(pair[0].id, pair[1].id, hd)
