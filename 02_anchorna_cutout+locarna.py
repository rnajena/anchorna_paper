# (C) 2025, Tom Eulenfeld, MIT license
from sugar import read, BioBasket
from anchorna import cutout, read_anchors
from glob import glob
import os
import time
from util import id_order_key


def create_input_fasta_files():
    seqs = read('comparison/pesti.fasta')
    anchors = read_anchors('comparison/anchors50.gff||8')
    N = len(anchors)
    assert N == 55
    anums = list(zip(range(-1, N), range(N+1)))
    for i in range(N+1):
        a1 = f'{i-1:02d}' if i > 0 else 'start'
        a2 = f'{i:02d}' if i < N else 'end'
        subseqs = cutout(seqs, anchors, a1 + '^' * (a1 != 'start'), a2 + '^' * (a2 != 'end'))
        if a1 == 'start':
            a1 = '000_start'
        subseqs.write(f'comparison/anchorna_cutout+locarna/{a1}_{a2}.fa')


def call_locarna_on_subsets():
    os.chdir('comparison/anchorna_cutout+locarna')
    t = time.time()
    for fname in sorted(glob('*.fa')):
        print(fname, flush=True)
        os.system(f'mlocarna -v --keep-sequence-order --stockholm --threads=2 {fname}')
    t2 = time.time()
    print('seconds used for LocARNA calls', t2-t)
    os.chdir('../..')


def merge_locarna_results():
    seqs = read('comparison/anchorna_cutout+locarna/*.out/results/result.stk')
    seqs.merge()
    seqs.str.replace('U', 'T')
    seqs.sort(id_order_key())
    seqs.write('comparison/pesti_anchorna_cutout+locarna.stk')


create_input_fasta_files()
call_locarna_on_subsets()
merge_locarna_results()
