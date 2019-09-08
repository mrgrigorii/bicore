import argparse
import os

from bicore import motif_tools
from bicore import utils


@utils.test_runtime
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--stats-folder',
        type=str,
        help='Folder with stats files',
    )
    parser.add_argument(
        '--output',
        type=str,
        help='output dir',
    )
    args = parser.parse_args()

    stats_files = os.listdir(args.stats_folder)

    stats = []

    for stats_file in stats_files:
        print(stats_file)
        with open(os.path.join(args.stats_folder, stats_file)) as fin:
            data = fin.read().strip('\n').split('\n')
            num_snps = len(data) - 1
            positions = int(data[0].split(':')[-1])
            name = data[0].split('|||')[0].lstrip('#')
            stats.append('\t'.join([
                name,
                str(num_snps),
                str(positions),
                str(round(num_snps / positions, 3)),
            ]))

    chromosome_data = motif_tools.load_fasta_from_gzip('chr1.fna.gz')
    positions = sum([1 for i in chromosome_data if i < 4])
    del chromosome_data[:]
    with open('sorted_snps_positions.txt', 'r') as fin:
        num_snps = len(fin.read().strip('\n').split('\n'))
        stats.append('\t'.join([
            'all_dna',
            str(num_snps),
            str(positions),
            str(round(num_snps / positions, 3)),
        ]))

    utils.dump_data('\n'.join(stats), os.path.join(args.output, 'ststs'))


if __name__ == '__main__':
    main()
