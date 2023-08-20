import sys
import logging
import math
from datetime import datetime, timedelta, date, time
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import namedtuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import spaceweather as sw

from ..stats.stats import Stats


RED = np.array([255, 204, 204]) / 255
ORANGE = np.array([255, 204, 154]) / 255
YELLOW = np.array([255, 255, 204]) / 255
GREEN = np.array([204, 255, 204]) / 255
BLUE = np.array([189, 227, 255]) / 255

GREEN2 = np.array([52, 253, 47]) / 255
YELLOW2 = np.array([254, 219, 49]) / 255
RED2 = np.array([252, 38, 28]) / 255


KP_COLORS_DISPLAY = {0: 'w',
                     1: 'w',
                     2: 'w',
                     3: 'w',
                     4: 'w',
                     5: BLUE,
                     6: GREEN,
                     7: YELLOW,
                     8: ORANGE,
                     9: RED}


KP_COLORS_DOCUMENT = {0: GREEN2,
                      1: GREEN2,
                      2: GREEN2,
                      3: GREEN2,
                      4: YELLOW2,
                      5: RED2,
                      6: RED2,
                      7: RED2,
                      8: RED2,
                      9: RED2}


KP_HATCH_DOCUMENT = {0: '',
                     1: '',
                     2: '',
                     3: '',
                     4: '/',
                     5: '\\',
                     6: '\\',
                     7: '\\',
                     8: '\\',
                     9: '\\'}


"""
Mapping between style identifier and the tuple of color and hatch
mappings.
"""
STYLE_MAP = {'display':  (KP_COLORS_DISPLAY,
                          False),
             'document': (KP_COLORS_DOCUMENT,
                          KP_HATCH_DOCUMENT)}


class IndexStats(namedtuple('IndiciesStats', 'kp dst')):
    pass


def plot_indices(date1,
                 date2,
                 fig=None,
                 style='display',
                 bar_lw=0.2,
                 dst_lw=1.75,
                 edgecolor=(0.1, 0.1, 0.1),
                 title_dt_format='%Y-%m-%d %H:%M',
                 stats=False):
    """
    """
    # gather Kp
    d1_kp = datetime.combine(date1, time(1, 30))
    d2_kp = datetime.combine(date2, time(0))

    kp = []
    kp_dt = []
    sw.update_data()
    df_3h = sw.ap_kp_3h()
    for dt in pd.date_range(d1_kp, d2_kp, freq='3H', inclusive='left'):
        kp.append(df_3h.loc[dt.strftime('%Y-%m-%d %H:%M:%S')].Kp)
        kp_dt.append(dt)
    d1_dst = datetime.combine(date1, time(0))
    d2_dst = datetime.combine(date2, time(0))
    years = range(d1_dst.year, d2_dst.year + 1)
    df = pd.concat([sw.omnie_hourly(x, cache=True) for x in years])
    I = (df.index >= d1_dst) & (df.index <= d2_dst)
    dst_I = df[I].Dst
    dst = dst_I.values
    dst_dt = dst_I.index
    # create plot
    N_days = (d2_dst - d1_dst).total_seconds() / 60 / 60 / 24
    if fig is None:
        fig = plt.figure(figsize=(11 * N_days / 6, 5))
    kp_colors, kp_hatch = STYLE_MAP[style]
    # Kp subplot
    ax1 = plt.subplot(111)
    left = np.array(kp_dt) - timedelta(hours=1.5)
    width = 3 / 24
    height = kp
    color = [kp_colors[int(math.floor(x))] for x in kp]
    hatch = [kp_hatch[int(math.floor(x))] if kp_hatch else None for x in kp]
    left = [x.to_pydatetime() for x in left]
    bars = plt.bar(left,
                   height,
                   width=width,
                   color=color,
                   linewidth=bar_lw,
                   edgecolor=[edgecolor] * len(left))
    for bar, hatch_i in zip(bars, hatch):
        bar.set_hatch(hatch_i)
    ax1.xaxis_date()
    plt.ylim(0, 9)
    plt.xlabel('UT')
    plt.ylabel('3-hour Kp index')
    # Dst subplot
    ax2 = ax1.twinx()
    plt.plot(dst_dt,
             dst,
             lw=dst_lw,
             marker='None',
             ls='-')
    plt.ylabel('Hourly DST [nT]')
    d1_str = datetime.strftime(dst_dt[0], title_dt_format)
    d2_str = datetime.strftime(dst_dt[-1], title_dt_format)
    title = 'GFZ $K_p$ and Dst Indices, {} to {}'.format(d1_str, d2_str)
    plt.title(title)
    if stats:
        index_stats = IndexStats(Stats(*kp),
                                 Stats(*dst))
        return index_stats, fig, ax1, ax2
    else:
        return fig, ax1, ax2


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Create plot of the Kp and Dst indices.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('pdf_fname',
                        type=str,
                        help='file to store plot')
    parser.add_argument('date1',
                        type=date.fromisoformat,
                        help='start date')
    parser.add_argument('date2',
                        type=date.fromisoformat,
                        help='end date')
    parser.add_argument('--style',
                        '-s',
                        type=str,
                        choices=sorted(STYLE_MAP),
                        default='display',
                        help='plot style (display is more colorful and meant for screen display whereas document is high contrast and uses hatches for both color and black-white interpretation and meant for use in publication)')
    args = parser.parse_args(argv[1:])

    plot_indices(args.date1,
                 args.date2,
                 style=args.style)

    plt.savefig(args.pdf_fname,
               bbox_inches='tight')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
