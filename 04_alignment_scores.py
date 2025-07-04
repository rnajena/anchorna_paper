# (C) 2025, Tom Eulenfeld, MIT license

from sugar import read, BioBasket
from sugar.data import submat
from statistics import mean, median


def gap_scores(seqs, addgap=-11, extentgap=-1):
    for seq in seqs:
        score = 0
        for match in seq.matchall('[-.]+'):
            n = len(match.group())
            score += addgap + (n - 1) * extentgap
        seq.meta.gapscore = score


def score_2seqs(seq1, seq2, sm):
    assert len(seq1) == len(seq2)
    return sum(
        sm[l1][l2]
        for l1, l2 in zip(seq1.data, seq2.data)
        if l1 not in '-.' and l2 not in '-.') + seq1.meta.gapscore + seq2.meta.gapscore


def score_alignment(seqs, sm='BLOSUM62', **kw):
    sm = submat(sm)
    gap_scores(seqs, **kw)
    scores = {(seq1.id, seq2.id): score_2seqs(seq1, seq2, sm)
              for i, seq1 in enumerate(seqs)
              for seq2 in seqs[i+1:]}
    return scores


PATH= 'comparison/'
FNAMES = [
    PATH + 'pesti_hand_curated.stk',
    PATH + 'pesti_locarna.stk',
    PATH + 'pesti_anchorna_export+locarna.stk',
    PATH + 'pesti_anchorna_cutout+locarna.stk',
    PATH + 'pesti_mafft.aln',
    PATH + 'pesti_anchorna+mafft.aln',
    PATH + 'pesti_dialign.aln',
    PATH + 'pesti_chaos+dialign.aln',
    PATH + 'pesti_anchorna+dialign.aln',
    #PATH + 'pesti_anchorna+dialign_allanchors.aln',
]

seqso = read(PATH + 'pesti.gff')
median_len = median(len(s) for s in seqso)

for i, fname in enumerate(FNAMES):
    seqs = read(fname)
    assert len(seqs) == 50
    if i > 0:  # manual changes in sequences in hand-curated alignment
        assert len(set([len(s) for s in seqs])) == 1
        for seqid in seqso.ids:
            # if len(seqso.d[seqid]) != len(seqs.d[seqid].copy().str.replace('-', '')):
            #     from IPython import embed
            #     embed()

            s1 = seqso.d[seqid]
            s2 = seqs.d[seqid].copy().str.replace('-', '').str.replace('.', '')
            for i in range(min(len(s1), len(s2))):
                if s1.data[i] != s2.data[i]:
                    print('warning', seqid, i, s1.data[i], s2.data[i], len(s1), len(s2))
                    BioBasket([s1, s2]).write(PATH + f'error_{seqid}.fasta')
                    break
            assert len(seqso.d[seqid]) == len(seqs.d[seqid].copy().str.replace('-', '').str.replace('.', ''))
            assert seqso.d[seqid].data == seqs.d[seqid].copy().str.replace('-', '').str.replace('.', '').data
    score_matrix = score_alignment(seqs, sm='NUC', addgap=-15, extentgap=-1)
    score = mean(score_matrix.values())
    nscore = score / median_len
    tool = fname.split('_')[-1].split('.')[0]
    rlen = round(len(seqs[0]), -2)
    print(f'{tool:>20} {nscore:.2f}/{rlen}')
    # print('           score: {:.2e}'.format(mean(score_matrix.values())))
    # print('normalized score: {:.2e}'.format(mean(score_matrix.values()) / median_len))
    # print('             len:', len(seqs[0]))
print('median len:', median_len)


# For reference, the NUC score:
"""NUC
#
# This matrix was created by Todd Lowe   12/10/92
#
# Uses ambiguous nucleotide codes, probabilities rounded to
#  nearest integer
#
# Lowest score = -4, Highest score = 5
#
    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2
T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2
G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2
C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2
S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1
W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1
R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1
Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1
K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1
M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1
B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1
V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1
H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1
D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1
N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1

"""
