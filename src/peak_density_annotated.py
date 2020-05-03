import argparse
import logging
import os
import gzip
import numpy as np
import subprocess
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from collections import defaultdict
from pathlib import Path


logger = logging.getLogger(__name__)

class MData:
    """
    Information for a heatmap plot. Each row has genomic coords attached.
    coords: 6-column matrix, akin to bedfile values.
    values: Values obtained from deeptools (by default this has 600 cols and
        as many regions as the input bedfile)
    limits: Visualization colors limits (1-98 nan-percentile of all the values)
    mean_vals: Column-wise means (for profile plots)
    mean_limits: Y-axis Limits for plotting.
    sum_value: Summary value for the whole matrix (right now col/row-wise mean)
    annotations: Matrix with the same number of rows as coords and values, where
        each column corresponds to a category annotation (yes/no format),
        where the order must match coords and values.

    xlabel: x axis label
    ylabel: y axis label
    """
    def __init__(self,
        coords,
        values,
        xlabel='',
        ylabel='',
        annotations=None,
        annotations_labels=None,
        remove_extreme=False):

        self.coords = coords
        self.values = values
        self.annotations = annotations

        self._validate_shapes()

        self.annotations_labels = annotations_labels

        self._validate_labels()

        # Using a mask in order to keep all the regions, so when another 
        # matrix is used to sort, the number of elements does not differ.
        self.mask = np.ones(self.values.shape[0], dtype=bool)

        if remove_extreme:
            self.remove_extreme_values()

        self.limits = self._compute_limits(self.values[self.mask, ])
        self.mean_vals = np.nanmean(self.values[self.mask, ], axis=0)
        self.mean_limits = [min(self.mean_vals), max(self.mean_vals)]

        self.sum_value = np.nanmean(self.mean_vals)
        self.xlabel = xlabel
        self.ylabel = ylabel


    def _validate_labels(self):
        if self.annotations_labels:
            if self.annotations.shape[1] != len(self.annotations_labels):
                msg = ('Number of annotations columns and number of labels '
                       'differ: {}, {} ({})'.format(self.annotations.shape[1],
                            len(self.annotations_labels),
                            str(self.annotations_labels))
                )
                logger.error(msg)
                raise ValueError(msg)


    def _validate_shapes(self):
        if self.coords.shape[0] != self.values.shape[0]:
            msg = 'Coords and values number of rows differ: {} {}'.format(
                self.coords.shape[0],
                self.values.shape[0]
            )
            logger.error(msg)
            raise ValueError(msg)

        if self.annotations is not None:
            if self.annotations.shape[0] != self.values.shape[0]:
                msg = ('Annotations must have the same number of rows as values'
                       ': {}, {}'.format(self.annotations.shape[0],
                                         self.values.shape[0])
                )
                logger.error(msg)
                raise ValueError(msg)


    def sort_values(self, order=None):
        if order is None:
            ref_avgs = np.nanmean(self.values, axis=1)
            order = ref_avgs.argsort()[::-1]

        self.coords = self.coords[order]
        self.values = self.values[order]
        self.mask = self.mask[order]
        if self.annotations is not None:
            self.annotations = self.annotations[order]

        return order

    def sort_by_mean(self, values):
        ref_avgs = np.nanmean(values, axis=1)
        my_avg = np.nanmean(self.values, axis=1)

        mean_val_list = np.array([(ref_avgs[i] + my_avg[i])/2.0 for i in range(len(ref_avgs))])
        order = mean_val_list.argsort()[::-1]

        self.coords = self.coords[order]
        self.values = self.values[order]
        self.mask = self.mask[order]
        if self.annotations is not None:
            self.annotations = self.annotations[order]


    def remove_extreme_values(self, percentile=99):
        rowmeans = np.nanmean(self.values, axis=1)
        zMax = np.nanpercentile(rowmeans, percentile)
        self.mask = (rowmeans[:] <= zMax)
        
    def nrows(self):
        return self.values.shape[0]


    @classmethod
    def from_gzmat(cls, gzfile, xlabel='', ylabel='', annot_file=None, remove_extreme=False):
        """
        Builds a MData object from a gzmat file (gz file produced by deeptools).
        """
        values, coords, header = cls._parse_matrix(gzfile,
            gzipped=True,
            header=False,
            skip=1)

        annotations = None
        labels = None

        if annot_file:
            annotations, labels = cls._parse_annotations(annot_file)

        return cls(coords, values, xlabel, ylabel,
            annotations=annotations,
            annotations_labels=labels,
            remove_extreme=remove_extreme)


    def validate_label_list(self, a_list):
        if self.annotations is not None:
            for a in a_list:
                if a not in self.annotations_labels:
                    msg = "{} label not in annotation list.".format(a)
                    logger.error(msg)
                    raise ValueError(msg)


    def filter_annotations(self, annot_label_list):
        if self.annotations is not None:
            self.validate_label_list(annot_label_list)
            index_list = [self.annotations_labels.index(a) for a in annot_label_list]
            new_annotations = self.annotations[:,index_list]
            return new_annotations

        return None


    @classmethod
    def _parse_annotations(cls, annotations):
        matrix, coords, header_labels = cls._parse_matrix(
            annotations,
            skip=0,
            header=True,
            gzipped=False,
            coords_col=4)

        # Annotations ONLY
        # matrix = [row[4:] for row in matrix[1:]]
        # return np.array(matrix, dtype=int), columns
        return matrix, header_labels


    @classmethod
    def _parse_matrix(cls, filename, sep='\t', skip=1, header=False, gzipped=True, coords_col=6):
        fi = None
        if gzipped:
            fi = gzip.open(filename)
        else:
            fi = open(filename)

        result = []
        line = fi.readline().rstrip()
        header_fields = None
        if header:
            header_fields = line.split(sep)[coords_col:]
            line = fi.readline().rstrip()

        for i in range(0, skip):
            line = fi.readline().rstrip()

        while line:
            fields = line.split()
            result.append(fields)
            line = fi.readline().rstrip()

        fi.close()
        values, coords =  cls._build_matrix(result, coords_col=coords_col)

        return  values, coords, header_fields


    @classmethod
    def _build_matrix(cls, mat, coords_col=6):
        result = [ m[coords_col:] for m in mat ]
        coords = [ m[:coords_col] for m in mat ]

        result_mat = np.array(result, dtype='float')
        coords_mat = np.array(coords)
        return result_mat, coords_mat


    def _compute_limits(self, matrix):
        zMin = None
        zMax = None

        matrix_flatten = None
        if zMin is None:
            matrix_flatten = matrix.flatten()
            # try to avoid outliers by using np.percentile
            zMin = np.nanpercentile(matrix_flatten, 1.0)
            if np.isnan(zMin):
                zMin = None

        if matrix_flatten is None:
            matrix_flatten = matrix.flatten()

        # try to avoid outliers by using np.percentile
        zMax = np.nanpercentile(matrix_flatten, 98.0)
        if np.isnan(zMax) or zMax < zMin:
            logger.warning("Max value is either smaller than min or zero: {}".format(zMax))
            zMax = None

        return zMin, zMax

    def get_split_annotation(self):
        """
        Splits a column in integer numbers. This was thought for chromhmm. It
        will return a list of labels called: annot_label.1, ... n and a matrix
        where each value is split to the corresponding column
        """
        old_mat = self.annotations[:,0]
        
        min_val = 1
        max_val = int(max(old_mat))

        # HACK: We only do this so far with chromhmm
        max_val = 15
        new_mat = np.zeros(shape=(len(old_mat), max_val), dtype=float)

        for i in range(1, max_val):
            new_mat[old_mat==i, i-1] = 1

        return new_mat, ['{}'.format(i) for i in range(min_val, max_val+1)]

    def set_zeros_to_nans(self, mat):
        mat[mat == 0] = np.nan
        return mat

    def plot(self, axes, fig, cmap, interpolation, flank, body,
            override_max=None, override_mean_max=None, override_mean_min=None,
            annotate_subset=None, signal_label='',
            loci_label='', title='', annot_title=''):
        
        long_label = 20
        very_long_label = 30

        cmap_plot=cmap
        hmax=self.limits[1]
        if override_max:
            # For global scale multiple plots
            hmax = override_max

        heatmap_ax = axes[1,0]
        heatmap = heatmap_ax.imshow(self.values[self.mask,],
            interpolation=interpolation,
            aspect='auto',
            origin='upper',
            vmax=hmax,
            cmap=cmap_plot)

        heatmap_ax.set_yticks([])
        ticks = [0.01, flank/10, (flank+body)/10, self.values.shape[1]-0.01]
        heatmap_ax.set_xticks(ticks)
        heatmap_ax.set_ylabel('{} ({})'.format(loci_label, self.values.shape[0]), labelpad=10)
       
        # heatmap_ax.set_xlabel("RPGC", fontsize=12)
        heatmap_ax.set_xticklabels(['-{}kb'.format(flank/1000),
            '',
            '',
            '+{}kb'.format(flank/1000)])

        cbar_ax = axes[2, 0]
        fig.colorbar(heatmap, cax=cbar_ax, alpha=1, ticks=[0, hmax], orientation="horizontal")

        cbar_ax.spines['top'].set_visible(False)
        cbar_ax.spines['left'].set_visible(False)
        cbar_ax.spines['right'].set_visible(False)
        cbar_ax.spines['bottom'].set_visible(False)
        cbar_ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        cbar_ax.xaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
        cbar_ax.set_xlabel("RPGC")

        # Mean plot
        # 10% of the margins for giving the plot space (legend fitting, hopefully)
        avg_limit_min = self.mean_limits[0]
        avg_limit_max = self.mean_limits[1] 
        if override_mean_max:
            avg_limit_max = override_mean_max

        if override_mean_min:
            avg_limit_min = override_mean_min

        margin_y = (avg_limit_max - avg_limit_min) / 10

        avgmin = avg_limit_min-margin_y
        avgmax = avg_limit_max+margin_y

        profile_ax = axes[0,0]
        intervals = range(0, len(self.mean_vals))
        profile_ax.plot(intervals, self.mean_vals, linewidth=1)
        profile_ax.set_ylim(avgmin, avgmax)
        profile_ax.set_title(signal_label)
        profile_ax.tick_params(axis='both', which='major')
        profile_ax.tick_params(axis='both', which='minor')
        profile_ax.spines['top'].set_visible(False)
        profile_ax.spines['right'].set_visible(False)
        profile_ax.spines['bottom'].set_visible(False)
        profile_ax.set_xticks([])

        # Annotation
        plt.xticks(rotation='vertical')
        annotations_ax = axes[1,1]
        if self.annotations:

            annotations_labels = self.annotations_labels
            annotations_labels = [label.replace('_', ' ') for label in annotations_labels]
            logger.debug("Annotating, finally: {}".format(','.join(annotations_labels)))

            to_plot = self.annotations[self.mask,]
            current_annotations_labels = annotations_labels
            if annotate_subset:
                to_plot = self.filter_annotations(annotate_subset)
        

            if len(annotations_labels) == 1 and max(to_plot) > 1:
                to_plot, labels = self.get_split_annotation()
                current_annotations_labels = labels

            vmax_annotations = np.amax(to_plot)
            not_all_zeros = np.any(to_plot)
            cmap_annot = 'Greys'
            if not_all_zeros:
                # This is done so zero does not count in the color scale
                to_plot = self.set_zeros_to_nans(to_plot)
                annot_map = annotations_ax.imshow(to_plot,
                                            interpolation='none',
                                            aspect='auto',
                                            origin='upper',
                                            cmap=cmap_annot,
                                            vmin=0,
                                            vmax=vmax_annotations,
                                            rasterized=True)


            annotations_ax.set_xticks([])
            annotations_ax.set_yticks([])
            annotations_ax.axis('on')
            annotations_ax.set_yticks([])
            annotations_ax.set_xticks(range(len(current_annotations_labels)))
            annotations_ax.set_xticklabels(current_annotations_labels, fontsize=8, rotation='vertical')
            annotations_ax.set_title(annot_title, fontsize=12)

        else:
            annotations_ax.axis('off')


def configure_plot(height, nannot):
    mpl.rcParams['axes.linewidth'] = 0.5
    mpl.rcParams['font.size'] = 18
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10

    cbar_width = 1
    single_heatmap_width = 3
    heatmap_width = single_heatmap_width
    
    # All annotations are plotted in one matrix
    annotations_width=min(single_heatmap_width, 2*nannot)

    if nannot > 0:
        heatmap_width += single_heatmap_width

    figwidth = heatmap_width

    cbar_height = 0.5
    average_plot_height = 2.5

    figheight = height + cbar_height + average_plot_height

    width_ratios = [single_heatmap_width, annotations_width]
    height_ratios = [average_plot_height, height, cbar_height]

    fig, axes_grid = plt.subplots(nrows=3, ncols=2, figsize=(figwidth, figheight),
        gridspec_kw={"width_ratios":width_ratios, "height_ratios":height_ratios})

    axes_grid[0,1].axis('off')
    axes_grid[2,1].axis('off')

    return fig, axes_grid

def parse_order(orderfile):
    with open(orderfile) as fi:
        result = [int(v)-1 for v in fi.readlines()]
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Plot a matrix precomputed by deeptools. Optionally,
            plot annotations together with it''',
        add_help=False)
        # epilog='.')

    parser.add_argument('--matrix',
        help='''Input matrix file (deeptools npz output).''',
        required=True)

    parser.add_argument('--sortby',
        help='''Use another matrix as reference to sort the matrix values.
            This is useful for side-by-side panels''', default=None)

    parser.add_argument('--height', help='Plot height', default=10, type=float)

    parser.add_argument('--annotations', help='Annotations to plot', default=None)
    parser.add_argument('--flank', help='Size of the flanking region (as in deeptools)', default=2500, type=float)
    parser.add_argument('--body', help='Size of the body region (as in deeptools)', default=0, type=float)
    parser.add_argument('--cmap', help='Heatmap colormap (as in used by matplotlib).', default="Reds")
    parser.add_argument('--title', help='Figure title', default='')
    parser.add_argument('--annot_title', help='Annotation title', default='Annotations')
    parser.add_argument('--xlabel', help='X axis label (signal)', default='signal')
    parser.add_argument('--ylabel', help='Y axis label (loci)', default='loci')
    parser.add_argument('--heatmap_max', help='Override color scale max, default is 98perc', default=None)
    parser.add_argument('--profile_limits', help='Override profile limits (example: 0.05:0.9)', default=None)
    parser.add_argument('--remove_top', help='Remove perc 99 (this is useful for profile plots)', action='store_true')
    parser.add_argument('--custom_order', help='''Override the sorting by a custom order that comes as text file.
        This is useful for plotting a matrix from a clustering made in another tool, for example. 
        Custom order must be a list of integers
        with length==matrix rows.''', default=None)
    parser.add_argument('--out', help='Output file', default='heatmap.pdf')


    # args: precomp matrix list + annotations
    args = parser.parse_args()

    matrix = args.matrix
    annot = args.annotations

    mdata = MData.from_gzmat(matrix, annot_file=annot, remove_extreme=args.remove_top) 
    sort_by = args.sortby

    if args.custom_order:
        order = parse_order(args.custom_order)
        mdata.sort_values(order)

    elif sort_by and sort_by != '-':
        ref_mat = MData.from_gzmat(sort_by)
        ref_order = ref_mat.sort_values()
        mdata.sort_values(order=ref_order)
    else:
        mdata.sort_values()    

    fig, ax = configure_plot(height=args.height, nannot=15)

    mean_max = None
    mean_min = None

    if args.profile_limits:
        mean_min, mean_max = args.profile_limits.split(':')
        mean_min = float(mean_min)
        mean_max = float(mean_max)

    heatmap_max = None
    if args.heatmap_max:
        heatmap_max = float(args.heatmap_max)

    mdata.plot(ax, fig, args.cmap, 'bilinear', args.flank, args.body,
        signal_label=args.xlabel, loci_label=args.ylabel, annot_title=args.annot_title,
        override_max=heatmap_max, override_mean_max=mean_max, override_mean_min=mean_min)

    fig.suptitle(args.title)
    fig.savefig(args.out, dpi=300)
    