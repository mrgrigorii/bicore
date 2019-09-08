import gzip
import os
import time


def test_runtime(func):
    start_time = time.time()

    def wrupped_func(*args, **kwargs):
        func(*args, **kwargs)
        runtime = round(time.time() - start_time, 2)
        print('Runtime: %s sec' % (str(runtime)))
    
    return wrupped_func


@test_runtime
def main():
    motif_files = os.listdir('dataD/')

    #sort_positions('snps_positions.txt')
    find_snp_in_motivs()
    


def write_only_snp_positions():
    with open('snps_positions.txt', 'w', encoding='utf-8') as out:
        N = 0
        for line in snp_iteritems('chr_1.txt.gz'):
            line = line.split('\t')
            #print(N, len(line), line[11])
            out.write(line[11] + '\n')
            N += 1
            if N >= 100000000:
                break
            elif N % 100000 == 0:
                print(str(N / 1000000) + ' %')


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


def find_snp_in_motivs():
    motif_data_list = []
    result_list = []
    motif_files = os.listdir('tests/motif_positions')
    print(motif_files)

    for motif in motif_files:
        motif_data = load_hocomoco_motif_data('tests/motif_positions/' + motif)
        list_pos = get_list_pos_motif(motif_data)
        motif_data_list.append([motif, list_pos])

        num_positions = get_motif_len(motif_data[0]) * len(list_pos)
        
        header = '#' + motif + '|||num_positions:' + str(num_positions)
        result_list.append([header])

    #GENE
    #list_pos = gen_pos_list

    with open('tests/sorted_snp_positions_test', 'r') as f:
        start = False
        data_lines = f.read().split('\n')
        len_data = len(data_lines)

        n = 0
        for line in data_lines:
            if n % 10000 == 0:
                print(str((100 * n) / len_data) + ' %')
            n += 1

            try:
                pos = int(line)
            except:
                print('Error ', line)
                continue

            pos_2 = pos #- 1   # my_pos = real_pos - 1
                
                
            for i_in_res, motif in enumerate(motif_data_list):

                remove_list = []
                for i in motif[1]:
                    if i[0]<=pos_2<=i[1]:
                        result_list[i_in_res].append(str(pos))
                    elif i[1] < pos:
                        if i[1] + 1000 < pos_2:
                            remove_list.append(i)
                        continue
                    else: 
                        break
                
                for i in remove_list:
                    motif[1].remove(i)


    for i, motif in enumerate(result_list):
        out = open(motif_data_list[i][0] + '.txt', 'w')
        out.write('\n'.join(motif))
        out.close()


def get_motif_len(motif_name):
    with open('tests/test_pwms_motifs.txt', 'r') as f:
        data = f.read().split('>')[1:]
        for motif in data:
            if motif_name in motif:
                #print(motif_name, motif)
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