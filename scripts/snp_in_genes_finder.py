import argparse
import os

from bicore import snp_finder
from bicore import utils


@utils.test_runtime
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--annotation',
        type=str,
        help='Gff gene annotation_file',
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
    )
    args = parser.parse_args()

    result = find_snp_in_genes_positions(args.annotation, args.snps_positions)
    utils.dump_data(
        result, os.path.join(args.output, 'snps_in_genes_forward_strand_v2')
    )


def find_snp_in_genes_positions(annotation, snp_positions):
    with open(annotation, 'r') as fin:
        gene_annotation = fin.read().split('\n')
        gene_positions = []
        num_positions = 0
        snp_positions = snp_finder.load_snp_positions(snp_positions)
        genes = 0
        for line in gene_annotation:
            if line.startswith('#') or not line:
                continue
            line = line.split()
            if line[0] == 'NC_000001.11' and line[1] == 'BestRefSeq' and \
                    line[2] == 'gene' and line[6] == '+':
                genes += 1
                start = int(line[3])
                end = int(line[4])
                num_positions += end - start
                gene_positions.append([start, end])
        print('Num genes %s' % genes)
        del gene_annotation[:]
        gene_positions.sort(key=lambda x: x[0])

        snp_positions = snp_finder.load_snp_positions(snp_positions)
        snp_gene_positions = snp_finder.find_snp_in_positions(
            gene_positions, snp_positions
        )

        header = '#genes|||num_positions:' + str(num_positions)

        result_list = [header]
        result_list.extend(snp_gene_positions)

        return '\n'.join(result_list)


if __name__ == "__main__":
    main()
