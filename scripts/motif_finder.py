import argparse
import os

from bicore import motif_tools
from bicore import utils


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fasta',
        type=str,
        help='Fasta file with one contig',
    )
    parser.add_argument(
        '--hocomoco',
        type=str,
        default='HOCOMOCOv11_core_pwms_HUMAN_mono.txt',
        help='hocomoco motifs file (pwms)',
    )
    parser.add_argument(
        '--output',
        type=str,
        help='output dir',
        default='data',
    )
    parser.add_argument(
        '--treshold',
        type=float,
        help='Motif quality treshold \% of max',
        default=0.9,
    )
    args = parser.parse_args()

    print("Start")
    motif_list = utils.load_hocomoco_motif(args.hocomoco)
    num_motifs = len(motif_list)
    fasta_data = motif_tools.load_fasta_from_gzip(args.fasta)

    print('Number of motifs: %s' % len(motif_list))
    num_motif = 1
    num_motifs = len(motif_list)

    for motif in motif_list:
        print('Run motif {}/{}'.format(num_motif, num_motifs))
        result = motif_tools.find_motif_in_fasta(
            fasta_data, motif, treshold=args.treshold
        )
        filename = os.path.join(args.output, '%s.hcm' % motif[0])
        utils.dump_data(result, filename)
        num_motif += 1


if __name__ == "__main__":
    main()
