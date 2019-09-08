import argparse
import collections
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
        '--motif-pos',
        type=str,
        help='Folder with motif positions files',
    )
    parser.add_argument(
        '--output',
        type=str,
        help='output dir',
    )
    args = parser.parse_args()

    stats_files = os.listdir(args.stats_folder)
    motiv_pos_files = os.listdir(args.motif_pos)

    # shuffled_motifs_pwms HOCOMOCOv11_core_pwms_HUMAN_mono.txt
    motifs_min_scores = get_motifs_min_score('shuffled_motifs_pwms')

    stats = []
    round_d = 1
    for stats_file in stats_files:
        print(stats_file)
        motiv_pos_file = None
        name = None
        motif_len = None
        stat = []
        min_score = get_min_score_for_motif(stats_file, motifs_min_scores)
        min_score_line = 'min_score ' + str(min_score)

        for motiv_pos in motiv_pos_files:
            # motiv_pos[:-13] delete tail AIRE_HUMAN.H11MO.0.C[.hcm.snpspos]
            if motiv_pos[:-13] in stats_file:
                motiv_pos_file = motiv_pos
                break
        if not motiv_pos_file:
            print('Error, not motiv_pos_file')
            continue

        motif_counter = collections.Counter()
        with open(os.path.join(args.motif_pos, motiv_pos_file)) as fin:
            data = fin.read().strip('\n').split('\n')
            motif_len = int(data[0].split(':')[1])
            data = [round(float(i.split(':')[1]), round_d) for i in data[1:]]
            for i in data:
                motif_counter[i] += 1

        snp_counter = collections.Counter()
        with open(os.path.join(args.stats_folder, stats_file)) as fin:
            data = fin.read().strip('\n').split('\n')
            name = data[0].split('|||')[0].lstrip('#')
            data = [round(float(i.split(':')[1]), round_d) for i in data[1:]]
            for i in data:
                snp_counter[i] += 1

        for key, value in snp_counter.items():
            positions = motif_len * motif_counter[key]
            num_snps = value
            if positions == 0:
                positions = 1

            stat.append([
                num_snps,
                positions,
                key,
                round(num_snps / positions, 2),
            ])

        stat.sort(key=lambda x: x[2])
        stat = ['\t'.join(list(map(str, i))) for i in stat]
        stat = '\n'.join(stat)
        print(stat)

        stats.append('{}\n{}\n{}'.format(name, min_score_line, stat))

    # shuffled_ststs_deep regular_ststs_deep
    utils.dump_data(
        '\n'.join(stats), os.path.join(args.output, 'all_с_shuffled_ststs_deep_075')
    )


def get_motifs_min_score(motifs_pwms_filename):
    motifs = utils.load_hocomoco_motif(motifs_pwms_filename)
    motifs_min_scores = {}
    for motif in motifs:
        min_score = utils.worst_motif_score(motif)
        motifs_min_scores[motif[0]] = min_score
    return motifs_min_scores


def get_min_score_for_motif(motif_stats_name, motifs_min_scores):
    for motif_base_name in motifs_min_scores.keys():
        if motif_base_name in motif_stats_name:
            return motifs_min_scores[motif_base_name]
    return 'Error, no motif in pwms'


if __name__ == '__main__':
    main()
