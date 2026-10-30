from itertools import combinations
import numpy as np


def _convert2np(seq):
    return np.frombuffer(seq.data.encode(), dtype='S1')


def spscore(ali, aliref):
    ids = set(aliref.ids)
    assert set(ali.ids) == ids
    m = 0
    n = 0

    for id1, id2 in combinations(sorted(ids), 2):
        sref1 = _convert2np(aliref.d[id1])
        sref2 = _convert2np(aliref.d[id2])
        s1 = _convert2np(ali.d[id1])
        s2 = _convert2np(ali.d[id2])
        assert len(sref2) == len(sref1)
        assert len(s2) == len(s1)
        gap_ref1 = (sref1 == b'-') | (sref1 == b'.')
        gap_ref2 = (sref2 == b'-') | (sref2 == b'.')
        gap_1 = (s1 == b'-') | (s1 == b'.')
        gap_2 = (s2 == b'-') | (s2 == b'.')

        # 1. Map every residue to its exact column index in the estimated alignment (ali)
        col_of_res1_in_ali = np.where(~gap_1)[0]
        col_of_res2_in_ali = np.where(~gap_2)[0]

        # 2. Find columns where both sequences have a residue in the reference (aliref)
        aligned_resids_ref = ~gap_ref1 & ~gap_ref2
        cols_ref = np.where(aligned_resids_ref)[0]

        if len(cols_ref) == 0:
            continue

        # 3. Find the 0-based residue indices for those reference columns
        r1_ref = np.cumsum(~gap_ref1)[cols_ref] - 1
        r2_ref = np.cumsum(~gap_ref2)[cols_ref] - 1

        # 4. Look up what columns these specific residues ended up in inside 'ali'
        c1_ali = col_of_res1_in_ali[r1_ref]
        c2_ali = col_of_res2_in_ali[r2_ref]

        # 5. If they share the same column in 'ali', the pair is successfully recovered
        m += np.sum(c1_ali == c2_ali)
        n += len(cols_ref)

    return m / n if n > 0 else 0.0

