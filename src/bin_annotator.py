import argparse
import os
import pandas as pd
import subprocess
import tempfile

def bedtools_intersect(filename_a, filename_b, out='annotations.bed', intersect_overlap=0):
    # bedtools intersect -wao -a "${peaks_file}" -b "${annotations_file}" > "${peaks_base}_${annot_base}.intersect.bed"
    intersect_tmp = tempfile.NamedTemporaryFile(prefix='intersect_', suffix='.bed').name
    groupby_tmp = tempfile.NamedTemporaryFile(prefix='groupby_', suffix='.bed').name
    annotations_tmp = tempfile.NamedTemporaryFile(prefix='annotations_', suffix='.bed', delete=False).name

    if intersect_overlap > 0:
        intersect_command = 'bedtools intersect -wao -a {} -b {} -F {} > {}'.format(
            filename_a,
            filename_b,
            str(intersect_overlap),
            intersect_tmp)
    else:
        intersect_command = 'bedtools intersect -wao -a {} -b {} > {}'.format(
            filename_a,
            filename_b,
            intersect_tmp)

    print(intersect_command)
    subprocess.check_call(intersect_command, shell=True)

    # bedtools groupby "${peaks_base}_${annot_base}.intersect.bed" -c 9,10 -o collapse,collapse > "${peaks_base}_${annot_base}.grouped.bed"
    # Bins file must have 3 columns. Annotations must have 6!
    groupby_command = 'bedtools groupby -i {} -c 7,10 -o collapse,collapse > {}'.format(
        intersect_tmp,
        annotations_tmp
    )
    print(groupby_command)
    subprocess.check_call(groupby_command, shell=True)

    select_max_overlap(annotations_tmp, out=out)
    os.remove(annotations_tmp)
    return out

def select_max_overlap(collapsed_file, id_col=3, val_col=4, sep='\t', out='annotations.bed'):
    fi = open(collapsed_file)
    fo = open(out, 'w')
    line = fi.readline().rstrip()

    within_field_sep = ','

    while line:
        fields = line.split(sep)
        ids = fields[id_col].split(within_field_sep)
        values = fields[val_col].split(within_field_sep)

        values_int = [int(i) for i in values]
        max_value = max(values_int)
        index = values_int.index(max_value)
        max_id = ids[index]

        new_line = sep.join(fields[0:id_col] + [max_id, values[index]] + fields[id_col+2:])
        fo.write(new_line + '\n')

        line = fi.readline().rstrip()

    fi.close()
    fo.close()

def merge_intersections(filelist):
    dflist = []
    for filename in filelist:
        df = pd.read_csv(filename, sep='\t', header=None)
        dflist.append(df)

    acum_df = dflist[0]
    for df in dflist[1:]:
        acum_df = acum_df.merge(df, on=[0,1,2])

    return acum_df

def main(args):
    # 1. Get only the bins coords and intersect them with the annotations.
    binfiles = args.binfiles.split(',')
    bedfiles = args.annotations.split(',')

    bins_tmp = tempfile.NamedTemporaryFile(
        prefix='bins_',
        suffix='.bed',
        delete=False).name

    # TODO: Check binfiles have the same binsize
    cut_command = 'cut -f 1-3 {} > {}'.format(binfiles[0], bins_tmp)
    subprocess.check_call(cut_command, shell=True)

    intersections = []
    for i, bed in enumerate(bedfiles):
        print("Intersecting {} with bins".format(bed))
        cur_file = bedtools_intersect(
            bins_tmp,
            bed,
            out='intersection_{}.bed'.format(i),
            intersect_overlap=float(args.minoverlap))

        print("Done")
        intersections.append(cur_file)

    annotations = merge_intersections(intersections)
    print("Annotations merged")
    # for i in intersections:
    #     os.remove(i)
    print(annotations)

    bin_dfs = merge_intersections(binfiles)
    if args.notscored:
        bin_dfs = bin_dfs[[0,1,2]]

    # 2. Join this result with the score files
    result_table = bin_dfs.merge(annotations, on=[0,1,2], how='inner')
    category_only = []
    print("Result table annotated")
    result_header = ['chr', 'start', 'end']
    if not args.notscored:
        for b in binfiles:
            result_header.append(os.path.basename(os.path.splitext(b)[0]))

        result_header.extend(['annotations', 'length'])
    else:
        for b in bedfiles:
            annot_name = os.path.basename(os.path.splitext(b)[0])
            category_only.append(annot_name)
            result_header.extend([annot_name, '{}_len'.format(annot_name)])

    result_table.columns = result_header

    result_header = ['bin_id'] + result_header
    result_table['bin_id'] = result_table.index
    result_table = result_table[result_header]
    result_table.to_csv(args.out, sep='\t', float_format='%.3f', header=True, index=False)

    # Filter only the categories
    if args.notscored:
        result_header = ['chr', 'start', 'end', 'bin_id'] + category_only
        result_table = result_table[result_header]
        result_table.to_csv(args.out, sep='\t', float_format='%.3f', header=True, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Take a list of binned files from deepTools output
            and a list of bedfiles to annotate them with. Generates
            a big table with all the columns from the binned files (scores)
            and annotations attached to it corresponding to the other bedfiles.
            Bins are also assigned an ID.''')

    parser.add_argument('--binfiles',
        help='Comma separated file list',
        required=True)

    parser.add_argument('--annotations',
        help='Comma separated bedfiles',
        required=True)

    parser.add_argument('--mode',
        help='''How to annotate:
        ALL: One column per category, yes/no.
        MAX: Only the feature that overlaps the most.''',
        default='ALL')

    parser.add_argument('--minoverlap',
        help='''Minimum percentage of overlap
            to be considered as a fraction of the annotated feature.''',
        default=0.8)

    parser.add_argument('--out', help='Output file', default='out.tab')
    parser.add_argument('--notscored',
        help=('Binfiles are bedfiles without scores, so only annotating a bed'
              'file with a list of bed files'),
        action='store_true')

    args = parser.parse_args()

    main(args)
