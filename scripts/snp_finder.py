import gzip
import os

from bicore import utils


@utils.test_runtime
def main():
    motif_files = os.listdir('dataD/')

    # sort_positions('snps_positions.txt')
    # find_snp_in_motivs()
    # find_snp_in_motivs_positions()
    # find_snp_in_genes_positions()


def write_only_snp_positions():
    with open('snps_positions.txt', 'w', encoding='utf-8') as out:
        n = 0
        for line in snp_iteritems('chr_1.txt.gz'):
            line = line.split('\t')
            # print(n, len(line), line[11])
            out.write(line[11] + '\n')
            n += 1
            if n >= 100000000:
                break
            elif n % 100000 == 0:
                print(str(n / 1000000) + ' %')


def snp_iteritems(file_path):
    with gzip.open(file_path, 'rb') as fin:
        data = [fin.readline().decode('utf-8') for i in range(7)]
        print('DATA header list', data)
        while data:
            data = fin.readline().decode('utf-8')
            yield data
        return None


def sort_positions(file_path):
    with open(file_path, 'r', encoding='utf-8') as fin:
        data = sorted([i for i in fin.read().split('\n') if i.strip(' ')])
        data = sorted(list(map(int, data)))
        data = list(map(str, data))
        with open('sorted_' + file_path, 'w', encoding='utf-8') as out:
            out.write('\n'.join(data).strip(' \n') + '\n')


def load_hocomoco_motif_data(file_name):
    with open(file_name, 'r') as f:
        data = f.read().strip('\n').split('\n')
        motif_name = data[0]
        data = [i.split(':') for i in data[1:]]
        data = [[int(i[0]), float(i[1])] for i in data]
    return motif_name, data


def find_snp_in_motivs_positions():
    motif_files = os.listdir('data')  #'tests/motif_positions')
    print(motif_files)
    snp_positions = load_snp_positions()

    for motif in motif_files:
        motif_data = load_hocomoco_motif_data('data/' + motif)
        motif_list_pos = get_list_pos_motif(motif_data)

        snp_motif_positions = find_snp_in_positions(
            motif_list_pos, snp_positions
        )

        num_positions = get_motif_len(motif_data[0]) * len(motif_list_pos)
        header = '#' + motif + '|||num_positions:' + str(num_positions)

        result_list = [header]
        result_list.extend(snp_motif_positions)
        print(len(result_list))

        with open('results/%s.txt' % motif, 'w') as out:
            out.write('\n'.join(result_list))
            out.close()


def find_snp_in_genes_positions():
    with open('GRCh38_latest_genomic.gff', 'r') as fin:
        gene_annotation = fin.read().split('\n')
        gene_positions = []
        num_positions = 0

        for line in gene_annotation:
            if line.startswith('#') or not line:
                continue
            line = line.split()
            if line[2] == 'gene':
                start = int(line[3])
                end = int(line[4])
                num_positions += end - start
                gene_positions.append([start, end])

        snp_positions = load_snp_positions()
        snp_gene_positions = find_snp_in_positions(
            gene_positions, snp_positions
        )

        header = '#genes|||num_positions:' + str(num_positions)

        result_list = [header]
        result_list.extend(snp_gene_positions)
        print(len(result_list))

        with open('results/genes.txt', 'w') as out:
            out.write('\n'.join(result_list))
            out.close()


def load_snp_positions():
    with open('sorted_snps_positions.txt', 'r') as f:
        data_lines = f.read().split('\n')
        data_lines = [int(i) for i in data_lines if i]
        return data_lines


def find_snp_in_positions(motif_list_pos, snp_positions):
    result_list = []
    len_data = len(snp_positions)

    n = 0
    start = 0
    for pos in snp_positions:
        if n % 1000000 == 0:
            print(str((100 * n) / len_data) + ' %')
        n += 1

        for i in range(start, len(motif_list_pos)):
            pos_vec = motif_list_pos[i]
            if pos_vec[1] < pos:
                start = i
                continue
            elif pos_vec[0] <= pos <= pos_vec[1]:
                result_list.append(str(pos))
            elif pos_vec[0] > pos:
                break

    return result_list


def get_motif_len(motif_name):
    with open('HOCOMOCOv11_core_pwms_HUMAN_mono.txt', 'r') as f:
        data = f.read().split('>')[1:]
        for motif in data:
            if motif_name in motif:
                # print(motif_name, motif)
                return len(motif.strip('\n').split('\n')) - 1


def get_list_pos_motif(motif_data, treshhold=0):
    motif_len = get_motif_len(motif_data[0])
    list_pos = []
    for i in motif_data[1]:
        if i[1] > treshhold:
            list_pos.append([i[0], i[0] + motif_len - 1])
    return list_pos 


if __name__ == '__main__':
    main()
