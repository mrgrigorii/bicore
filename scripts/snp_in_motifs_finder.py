import argparse
import os

from bicore import snp_finder
from bicore import utils


@utils.test_runtime
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--motifs',
        type=str,
        help='Folder that contain motif positions files',
        default='updated_results'
    )
    parser.add_argument(
        '--snps-positions',
        type=str,
        help='Sorted snps positions file',
    )
    parser.add_argument(
        '--output',
        type=str,
        help='output dir',
    )
    args = parser.parse_args()

    snp_positions = snp_finder.load_snp_positions(args.snps_positions)
    motif_files = os.listdir(args.motifs)
    print(motif_files)

    for motif in motif_files:
        print('Motif: ' + str(motif))
        result = find_snp_in_motiv_positions(args.motifs, motif, snp_positions)
        filename = os.path.join(args.output, '%s.snpspos' % motif)
        dump_large_data(result, filename)


def dump_large_data(data, filename):
    utils.check_folder_exists(filename)
    with open(filename, 'w') as out:
        len_result = len(data)
        part_len = len_result // 3
        print('Dumping')
        out.write('\n'.join(data[0:part_len]))
        out.write('\n')
        out.write('\n'.join(data[part_len:2 * part_len]))
        out.write('\n')
        out.write('\n'.join(data[2 * part_len:]))


def load_hocomoco_motif_data(file_name):
    with open(file_name, 'r') as f:
        data = f.read().strip('\n').split('\n')
        motif_name, motif_len = data[0].split(':')
        motif_len = int(motif_len)
        data = [i.split(':') for i in data[1:]]
        data = [[int(i[0]), float(i[1])] for i in data]
    return motif_name, motif_len, data


def find_snp_in_motiv_positions(motifs_folder, motif, snp_positions):
    motif_data = load_hocomoco_motif_data(
        os.path.join(motifs_folder, motif)
    )
    motif_list_pos = get_list_pos_motif(motif_data)

    snp_motif_positions = snp_finder.find_snp_in_positions(
        motif_list_pos, snp_positions, is_motif=True,
    )

    num_positions = motif_data[1] * len(motif_list_pos)
    header = '#' + motif + '|||num_positions:' + str(num_positions)

    result_list = [header]
    result_list.extend(snp_motif_positions)
    return result_list  # '\n'.join(result_list)


def get_list_pos_motif(motif_data, treshhold=0):
    motif_len = motif_data[1]
    list_pos = []
    for i in motif_data[2]:
        if i[1] > treshhold:
            list_pos.append([i[0], i[0] + motif_len - 1, i[1]])
    return list_pos


if __name__ == '__main__':
    main()
