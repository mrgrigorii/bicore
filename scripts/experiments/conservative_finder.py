import argparse
import collections
import gzip
import os

from bicore import snp_finder
from bicore import utils


@utils.test_runtime
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fasta',
        type=str,
        help='Fasta file with one contig',
    )
    parser.add_argument(
        '--snps-positions',
        type=str,
        help='Sorted snps positions file',
        default='sorted_snps_positions.txt',
    )
    parser.add_argument(
        '--output',
        type=str,
        help='output dir',
        default='data',
    )
    args = parser.parse_args()

    print("Start")
    fasta_data = load_fasta_from_gzip(args.fasta)
    snp_positions = snp_finder.load_snp_positions(args.snps_positions)

    motifs = collections.Counter()

    len_data = len(fasta_data)
    len_motif = 8
    current_snp = 0

    for i in range(0, len_data - len_motif):
        if i % 5000000 == 0:
            progress = str(round(i * 100.0 / len_data, 2))
            print(progress + ' %')

        for j in range(current_snp, len(snp_positions)):
            if snp_positions[j] < i and j > current_snp:
                current_snp = j
            elif i <= snp_positions[j] <= i + len_motif:
                motifs[fasta_data[i:i + len_motif]] += 1
            elif snp_positions[j] > i + len_motif:
                break

    motifs_list = list(motifs.items())
    motifs_list.sort(key=lambda x: x[1])
    motifs_list = ['\t'.join(list(map(str, i))) for i in motifs_list]
    utils.dump_data(
        '\n'.join(motifs_list),
        os.path.join(args.output, 'cons_motifs_stat')
    )


def load_fasta_from_gzip(fasta_gzip_name):
    with gzip.open(fasta_gzip_name, 'rb') as f:
        data = f.read().decode("utf-8")
        data = ''.join(data.strip('\n').split('\n')[1:])
        return data


if __name__ == "__main__":
    main()
