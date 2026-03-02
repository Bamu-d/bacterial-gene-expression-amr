#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd

from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Patch

def get_options():
    description = "Plot a sequence 'alignment' a-la panfeed"

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("input",
            help="Input FASTA file")
    parser.add_argument("output",
            help="Output name, without the file extension")

    parser.add_argument('--random-sample',
                        type=float,
                        default=None,
                        help="Only show a randomly picked set of sample "
                             "(a value between 0 and 1 indicating the proportion"
                             " to show, default all)")
    parser.add_argument('--samples',
                        default=None,
                        help="File with one sample to show per line, "
                             "(the same order as in the file will be shown,"
                             " default all)")

    parser.add_argument('--upstream',
                        type=int,
                        default=200,
                        help="Number of bases upstream of gene start "
                             "(default 200)")
    parser.add_argument('--downstream',
                        type=int,
                        default=30,
                        help="Number of bases downstream of gene start "
                             "(default 30)")
    parser.add_argument('--start',
                        type=int,
                        default=None,
                        help="Relative position to start the plots "
                             "(default all available positions)")
    parser.add_argument('--stop',
                        type=int,
                        default=None,
                        help="Relative position to end the plots "
                             "(default all available positions)")

    parser.add_argument('--format',
                        choices=('png',
                                 'tiff',
                                 'pdf',
                                 'svg'),
                        default='png',
                        help='Output format for plots (default %(default)s)')

    parser.add_argument('--dpi',
                        type=int,
                        default=300,
                        help='Output resolution (DPI, default %(default)d)')

    parser.add_argument('--xticks',
                        type=int,
                        default=50,
                        help="Spacing for ticks on x axis (default %(default)d)")
    parser.add_argument('--print-samples',
                        action='store_true',
                        default=False,
                        help="Write sample names on the y-axis (default: don't)")

    parser.add_argument('--height',
                        type=float,
                        default=9.,
                        help="Figure height (inches, default %(default).1f)")

    parser.add_argument('--width',
                        type=float,
                        default=10.,
                        help="Figure width (inches, default %(default).1f)")

    parser.add_argument('--hide-bases',
                        action='store_true',
                        default=False,
                        help="Hide the individual base letters in the plot (default: show them)")

    return parser.parse_args()


def main():
    args = get_options()

    if args.random_sample is not None and args.samples is not None:
        sys.stderr.write("only one argument between --samples and --random-sample should be used\n")
        sys.exit(1)

    if args.samples is not None:
        samples = [x.rstrip() for x in open(args.samples)]

    r = []
    s = ''
    sid = ''

    base2int = {'A': 0,
                'G': 1,
                'T': 2,
                'C': 3}

    for l in open(args.input):
        if l.startswith('>') and s != '':
            for i, c in enumerate(s.upper()):
                if c not in 'AGCT':
                    c = 'N'
                r.append((sid, i, c, base2int.get(c, np.nan)))
            s = ''
            sid = l.rstrip()[1:].split(':')[0].split('_')[0]
        elif l.startswith('>'):
            sid = l.rstrip()[1:].split(':')[0].split('_')[0]
            continue
        else:
            s += l.rstrip()
    for i, c in enumerate(s.upper()):
        if c not in 'AGCT':
            c = 'N'
        r.append((sid, i, c, base2int.get(c, np.nan)))

    r = pd.DataFrame(r, columns=['sample', 'start', 'base', 'scalar'])
    
    # Convert positions to relative coordinates
    r['start'] = r['start'] - args.upstream
    
    if args.start is not None:
        r = r[(r['start'] >= args.start) &
              (r['start'] <= args.stop)]

    def handle_paralogs(x):
        if len(x) > 1:
            return 99
        return x

    def handle_paralogs_text(x):
        if len(x) > 1:
            return '-'
        return x

    base_colors = list(sns.color_palette('tab20', 4))
    cmap = colors.LinearSegmentedColormap.from_list('nucleotides', base_colors, 4)
    cmap.set_bad('xkcd:grey')
    cmap.set_over(list(sns.color_palette('tab20', 5))[-1])

    legend_elements = [
        Patch(facecolor=base_colors[0], edgecolor='black', label='A'),
        Patch(facecolor=base_colors[1], edgecolor='black', label='G'),
        Patch(facecolor=base_colors[2], edgecolor='black', label='T'),
        Patch(facecolor=base_colors[3], edgecolor='black', label='C'),
        Patch(facecolor='grey', edgecolor='black', label='N')
    ]

    b = r.pivot_table(index='sample', columns='start',
                      values='scalar',
                      aggfunc=handle_paralogs)
    t = r.pivot_table(index='sample', columns='start',
                      values='base',
                      aggfunc=handle_paralogs_text)

    if args.samples is not None:
        b = b.loc[samples]
        t = t.loc[b.index]
    if args.random_sample is not None:
        b = b.sample(frac=args.random_sample)
        t = t.loc[b.index]

    # Create x-axis ticks at specified intervals
    start_pos = -args.upstream
    end_pos = args.downstream
    xtick_positions = np.arange(start_pos, end_pos + 1, args.xticks)
    
    # Create mapping between matrix positions and relative positions
    pos_mapping = {i: pos for i, pos in enumerate(range(start_pos, end_pos + 1))}
    
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    im = ax.imshow(b.values, cmap=cmap, vmin=0, vmax=3,
                   aspect="auto", interpolation='none',
                   rasterized=True, alpha=1)
    
    if not args.hide_bases:
        for x in range(t.shape[1]):
            for y in range(t.shape[0]):
                ax.text(x, y, t.iloc[y, x],
                        ha='center', va='center')

    if args.print_samples:
        ax.set_yticks(range(b.shape[0]), b.index)
    else:
        ax.set_yticks([])

    # Set x-axis ticks to show relative positions
    ax.set_xticks(np.where(np.isin(list(pos_mapping.values()), xtick_positions))[0])
    ax.set_xticklabels(xtick_positions)
    
    # Add vertical line at gene start (position 0)
    gene_start_idx = np.where(np.array(list(pos_mapping.values())) == 0)[0][0]
    ax.axvline(gene_start_idx - 0.5, lw=3, color="black")
    
    ax.set_xlabel("position relative to gene start")
    ax.set_ylabel("samples")

    ax.legend(handles=legend_elements, loc='upper right', title="Nucleotides")

    plt.savefig(f'{args.output}.{args.format}', dpi=args.dpi, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()