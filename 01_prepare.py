# (C) 2025, Tom Eulenfeld, MIT license
from sugar import read
from util import get_50_ids


def prepare_50_seqs():
    seqs = read('comparison/pesti_example.gff')
    seqs.select(id_in=get_50_ids(), inplace=True)
    seqs.write('comparison/pesti.gff')
    seqs.write('comparison/pesti.fasta')


prepare_50_seqs()
