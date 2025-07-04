# (C) 2025, Tom Eulenfeld, MIT license
from sugar import read
from util import get_50_ids, id_order_key


def postprocess_locarna():
    seqs = read('comparison/pesti_hand_curated_55seqs_raw.stk')
    for seq in seqs:
        seq.id = seq.id.split('_', 2)[-1]
    seqs.sort(id_order_key())
    seqs.str.replace('U', 'T')
    seqs.write('comparison/pesti_hand_curated_55seqs.stk')
    seqs.select(id_in=get_50_ids(), inplace=True)
    seqs.sort(id_order_key())
    seqs.write('comparison/pesti_hand_curated.stk')

    seqs = read('comparison/pesti_locarna_raw.stk')
    seqs.str.replace('U', 'T')
    seqs.sort(id_order_key())
    seqs.write('comparison/pesti_locarna.stk')

    seqs = read('comparison/pesti_anchorna_export+locarna_raw.stk')
    seqs.str.replace('U', 'T')
    seqs.sort(id_order_key())
    seqs.write('comparison/pesti_anchorna_export+locarna.stk')


# now done on command line
# def merge_mafft_alis():
#     seqs1 = read('comparison/pesti4mafft_start_A14.aln')
#     seqs2 = read('comparison/pesti4mafft_A14_end.aln')
#     seqs = (seqs1 + seqs2).merge()
#     seqs.write('comparison/pesti_anchorna+mafft.aln', 'fasta')


def correct_dialign():
    fexpr = 'comparison/pesti_{}.aln'
    for job in ('dialign', 'chaos+dialign', 'anchorna+dialign'):
        seqs = read(fexpr.format(job + '_uncorrected'), encoding='ISO-8859-15')
        seqs.str.replace('!0', '--')  # seqid KP343640
        seqs.str.replace('±.', '--')  # seqid C_077026
        seqs.write(fexpr.format(job), 'fasta')


postprocess_locarna()
correct_dialign()
