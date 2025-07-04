# (C) 2025, Tom Eulenfeld, MIT license

from anchorna import read_anchors, cutout
from collections import defaultdict
import matplotlib.pyplot as plt
from sugar import read, read_fts
PATH = 'comparison/'
FIGPATH = 'figs/'


def _get_arrow_pos(cseqs):
    a3pos = set()
    for s in cseqs:
        start, stop, _ = s.slindex(gap='-.')[s.fts[3]].indices(len(s))
        a3pos.add((start + stop) // 2)
    return a3pos


def _anchors_groupby(ft):
    return ft.name.split('_')[0]


def full_alignment_plot():
    seqs = read(PATH + 'pesti_hand_curated_55seqs.stk')
    anchors = read_anchors(PATH + 'anchors.gff||8')
    seqs.fts = anchors.convert2fts(mode='nt')
    for ft in seqs[-1].fts:  # Sandra cut of the 5'UTR
        ft.loc.start -= 518
        ft.loc.stop -= 518
        #f7c743
    ax = seqs.plot_alignment(fts=True, fts_color='C1', fts_alpha=0.8, color='0.4', figsize=(16, 3), dpi=600, rasterized=True)
    ax.get_figure().savefig(FIGPATH + 'full_alignment_hand_curated.png')
    seqs = read(PATH + 'pesti_locarna.stk')
    seqs.plot_alignment(FIGPATH + 'full_alignment_locarna.png', figsize=(16, 3), dpi=600)
    seqs = read(PATH + 'pesti_anchorna_cutout+locarna.stk')
    seqs.plot_alignment(FIGPATH + 'full_alignment_anchorna_cutout+locarna.png', figsize=(16, 3), dpi=600)
    seqs = read(PATH + 'pesti_anchorna_export+locarna.stk')
    seqs.plot_alignment(FIGPATH + 'full_alignment_anchorna_export+locarna.png', figsize=(16, 3), dpi=600)



def subalignment_plot():
    fnames = [
        PATH + 'pesti_hand_curated.stk',
        PATH + 'pesti_locarna.stk',
        PATH + 'pesti_anchorna_cutout+locarna.stk',
        PATH + 'pesti_mafft.aln',
        PATH + 'pesti_anchorna+mafft.aln',
        PATH + 'pesti_dialign.aln',
        PATH + 'pesti_chaos+dialign.aln',
        PATH + 'pesti_anchorna+dialign.aln',
    ]
    labels = [
        'hand-curated',
        'vanilla LocARNA',
        'AnchoRNA+LocARNA',
        'vanilla MAFFT',
        'AnchoRNA+MAFFT',
        'vanilla Dialign',
        'CHAOS+Dialign',
        'AnchoRNA+Dialign',
        ]

    anchors_orig = read_anchors(PATH +'anchors50.gff|11:18')
    anchors_handcurated = anchors_orig.copy()
    # The sequence OM030319 is altered in the hand-curated alignment.
    # Its 5' UTR was removed.
    # Corresponding anchors need to be adapted
    for anchor in anchors_handcurated:
        for fluke in anchor:
            if fluke.seqid == 'OM030319':
                fluke.offset = 0

    fig, axes = plt.subplots(len(fnames), figsize=(10, 10))
    for ax, label, fname in zip(axes, labels, fnames):
        anchors = anchors_handcurated if label == 'hand-curated' else anchors_orig
        anchors_nt = anchors.convert2fts(mode='nt')
        seqs = read(fname)
        # from util import id_order_key
        # seqs.sort(id_order_key())
        seqs.fts = anchors_nt
        cseqs = cutout(seqs, anchors, '0', '6', gap='-.', update_fts=True)
        cseqs.plot_alignment(ax=ax, color='0.7',
                             fts=True, fts_colorby=_anchors_groupby, fts_color_gap_alpha=0.3,
                             xticks=False, rasterized=True)
        ax.annotate(label, (0, 1.02), xycoords='axes fraction', size='small')
        ax.annotate(f'{round(len(cseqs[0]), -1)}nt', (0.3, 1.02), xycoords='axes fraction', size='small')
        for a3p in _get_arrow_pos(cseqs):
            ax.annotate('', (a3p, len(cseqs)), (a3p, len(cseqs)+8),
                        arrowprops=dict(arrowstyle='->', color = 'C3'), annotation_clip=False)
            if  label == 'AnchoRNA+Dialign':  # only display one arrow here
                break
    kw = dict(ha='center', xycoords='axes fraction', size='small')
    axes[0].annotate('nucleotides\n(gray)', (0.765, 1.02), **kw)
    axes[0].annotate('gaps\n(white)', (0.848, 1.02), **kw)
    axes[0].annotate('AnchoRNA anchors\n(colored)', (0.96, 1.02), **kw)
    fig.savefig(FIGPATH + 'alis.pdf', dpi=300)


def four_anchors_plot():
    from util import get_24_ids
    anchors = read_anchors(PATH +'anchors.gff|11:18')
    seqs = read(PATH + 'pesti.gff')
    seqs.fts = anchors.convert2fts(mode='nt')
    seqs = seqs.select(id_in=get_24_ids())
    cseqs = cutout(seqs, anchors, '0-9', '2+9', gap='-.', update_fts=True)
    cseqs2 = cutout(seqs, anchors, '3-18', '3+18', gap='-.', update_fts=True)
    cseqs = (cseqs + cseqs2).merge(spacer=9 * ' ', update_fts=True)
    print('removed nts between anchors:', cutout(seqs, anchors, '2>+9', '3<-18', gap='-.', update_fts=True).str.replace('-', '').tostr(h=None))

    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111)
    kw = dict(ax=ax, xticks=False, fts_display='box', fts_color=None, fts=True, fts_alpha=0.6,
            fts_colorby=_anchors_groupby, fts_box_groups=_anchors_groupby,
            symbols=True, symbol_size=6, edgecolors='white', linewidth=1, alpha=0.4)
    cseqs.plot_alignment(extent=[0, 50, -20, -2], color=None, **kw)
    tseqs = cseqs.copy().translate(complete=True, update_fts=True)
    tseqs.str.replace('X', ' ')
    tseqs.plot_alignment(extent=[0, 50, 0, 18], color='flower', **kw)
    ax.set_ylim(-20.5, 18.5)
    akw = dict(size='xx-large', annotation_clip=False)
    ax.annotate('A', (-0.8, 17.2), **akw)
    ax.annotate('B', (-0.8, -2.8), **akw)
    fig.savefig(FIGPATH + 'anchor_example.pdf', bbox_inches='tight')


def genome_anchors_plot():
    id_ = 'NC_076029'
    fts = read_anchors(PATH + 'anchors.gff').convert2fts(mode='nt')
    fts = fts.select(seqid_eq=id_)
    seqs = read(PATH + 'pesti.gff').select(id_eq=id_)
    seqs.fts = fts
    for ft in seqs.fts:
        if len(ft) < 30:
            ft.loc.stop = ft.loc.start + 30
    fts_color = defaultdict(lambda: (0, 0, 0))
    for i in range(7):
        fts_color[f'A1{i+1}'] = f'C{i}'
    seqs.plot_alignment(fts=True, color='white', alpha=0, fts_alpha=1, fts_colorby=_anchors_groupby, fts_color=fts_color,
                        xticks=False, figsize=(5,0.3), fname=FIGPATH + 'genome_anchors.png', transparent=True)


def genome_anchor_distances():
    anchors = read_anchors(PATH + 'anchors.gff')
    a0 = anchors[0].todict_seqid()
    min_ = 100
    max_ = 0
    for i, a in enumerate(anchors[1:]):
        a = a.todict_seqid()
        difs = [(a[id_].start - a0[id_].stop) for id_ in a0]
        a0 = a
        print(f'min/max distance between anchor {i:02d} and {i+1:02d}:', min(difs), max(difs))
        min_ = min(difs + [min_])
        max_ = max(difs + [max_])
    print('min distance between two anchors:', min_)
    print('max distance between two anchors:', max_)

    a13 = anchors[13].todict_seqid()
    a15 = anchors[15].todict_seqid()
    difs = [(a15[id_].start - a13[id_].stop) for id_ in a13]
    print('min distance between anchors A13 and A15:', min(difs))
    print('max distance between anchors A13 and A15:', max(difs))


subalignment_plot()
full_alignment_plot()
four_anchors_plot()
genome_anchors_plot()
genome_anchor_distances()
