import array
import os
import gzip


NUCL_STRING = "ACGT"


def load_fasta_from_file(fasta_file_name):
    with open(fasta_file_name, 'r') as f:
        data = f.read()
        return load_fasta_data(data)

def load_fasta_from_gzip(fasta_gzip_name):
    with gzip.open(fasta_gzip_name, 'rb') as f:
        data = f.read().decode("utf-8")
        return load_fasta_data(data)


def load_fasta_data(data):
    prep_data = array.array('B')
    data = ''.join(data.strip('\n').split('\n')[1:])
    for i in data:
        if i == NUCL_STRING[0]: prep_data.append(0)
        elif i == NUCL_STRING[1]: prep_data.append(1)
        elif i == NUCL_STRING[2]: prep_data.append(2)
        elif i == NUCL_STRING[3]: prep_data.append(3)
        else: prep_data.append(4)
    # should be Data len:  248956422
    print('Data len: ', len(prep_data))
    return prep_data


def load_hocomoco_motif(file_name):
    motif_list = []
    with open(file_name, 'r') as f:
        data = f.read().split('>')[1:]
        for motif in data:
            motif = motif.split('\n')
            motif_list.append([motif[0], list(list(map(float, i.split())) + [0.0] for i in motif[1:] if i)])

    return motif_list


def best_motif_score(motif):
    best_score = 0
    for i in motif[1]:
        best_score += max(i)
    return best_score


def find_motif_in_fasta(data, motif_list, treshhold=0, data_folder='/'):
    """ result format

    >motif_name
    position_0:treshold_0
    position_1:treshold_1
    ...
    """

    print(len(motif_list))
    num_motif = 0
    num_motifs = len(motif_list)

    for motif in motif_list:
        max_score = best_motif_score(motif)
        result = []
        num_motif += 1

        result.append(motif[0])
        len_motif = len(motif[1])

        len_data = len(data)

        for i in range(0, len_data - len_motif):
            if i%5000000 == 0:
                print(str(round(i*100.0/len_data, 2)) + '% motif ' + str(num_motif) + ' of ' + str(num_motifs))
                print(len(result))

            motif_trashhold = 0.0
            for j in range(0, len_motif):
                motif_trashhold += motif[1][j][data[i+j]]

            if motif_trashhold > max_score*0.5:
                result.append(str(i) + ':' + str(motif_trashhold))


        with open(data_folder + 'NEWchr2_' + result[0] + '.hcm', 'w') as out:
            out.write('\n'.join(result))


def main():
    print("Start")
    file_name = "HOCOMOCOv11_core_pwms_HUMAN_mono.txt"
    motif_list = load_hocomoco_motif(file_name)
    #print(best_motif_score(motif_list[0]))

    #find_motif_in_fasta(fasta_data, motif_list, treshhold=10**17)
    #data_folder = "/home/mrgrigorii/Projects/GRCh38.p12"
    fasta_data = load_fasta_from_gzip('chr1.fna.gz')
    result = find_motif_in_fasta(fasta_data, motif_list[110:], treshhold=0, data_folder='data/')


if __name__ == "__main__":
    main()
        
