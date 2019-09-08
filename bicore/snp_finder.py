import gzip


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


def load_snp_positions(snp_positions):
    with open(snp_positions, 'r') as f:
        data_lines = f.read().split('\n')
        data_lines = [int(i) for i in data_lines if i]
        return data_lines


def find_snp_in_positions(motif_list_pos, snp_positions, is_motif=False):
    result_list = []
    len_data = len(snp_positions)

    n = 0
    start = 0
    for pos in snp_positions:
        if n % 1000000 == 0:
            print(str(round((100 * n) / len_data, 2)) + ' %')
        n += 1

        for i in range(start, len(motif_list_pos)):
            pos_vec = motif_list_pos[i]
            if pos_vec[1] < pos:
                start = i
                continue
            elif pos_vec[0] <= pos <= pos_vec[1]:
                if is_motif:
                    result_list.append(
                        '{}:{}'.format(pos, pos_vec[2])
                    )
                else:
                    result_list.append(str(pos))
            elif pos_vec[0] > pos:
                break

    return result_list
