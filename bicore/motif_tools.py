import array
import gzip

from bicore import const
from bicore import utils


def load_fasta_from_file(fasta_file_name):
    with open(fasta_file_name, 'r') as f:
        data = f.read()
        return load_fasta_data_as_digital_array(data)


def load_fasta_from_gzip(fasta_gzip_name):
    with gzip.open(fasta_gzip_name, 'rb') as f:
        data = f.read().decode("utf-8")
        return load_fasta_data_as_digital_array(data)


def load_fasta_data_as_digital_array(data):
    ''' Convert ACGT to 0123'''
    prep_data = array.array('B')
    data = ''.join(data.strip('\n').split('\n')[1:])
    for i in data:
        if i == const.NUCL_STRING[0]:
            prep_data.append(0)
        elif i == const.NUCL_STRING[1]:
            prep_data.append(1)
        elif i == const.NUCL_STRING[2]:
            prep_data.append(2)
        elif i == const.NUCL_STRING[3]:
            prep_data.append(3)
        else:
            prep_data.append(4)

    return prep_data


def find_motif_in_fasta(data, motif, treshold=0.9):
    ''' result format

    >motif_name:len_motif
    position_0:treshold_0
    position_1:treshold_1
    ...
    '''
    max_score = utils.best_motif_score(motif)
    min_score = utils.worst_motif_score(motif)
    treshold_score = (max_score - min_score) * treshold + min_score
    print('Treshhold_score', treshold_score, min_score, max_score)
    len_motif = len(motif[1])
    len_data = len(data)

    result = []
    header = '#{name}:{length}\n'.format(
        name=motif[0],
        length=len_motif,
    )

    for i in range(0, len_data - len_motif):
        if i % 5000000 == 0:
            progress = str(round(i * 100.0 / len_data, 2))
            print(progress + ' %')
            print('Find %s motifs' % len(result))

        motif_score = 0.0
        for j in range(0, len_motif):
            motif_score += motif[1][j][data[i + j]]

        if motif_score > treshold_score:
            result.append([i, round(motif_score, 3)])

    filtered_result = motif_filter(result, len_motif)

    return header + filtered_result


def motif_filter(motif_positions, len_motif):
    '''Choose motifs with the best score from overlapping'''
    if not motif_positions:
        return ''
    filtered_positions = []
    current_pos = motif_positions[0]
    for pos in motif_positions[1:]:
        if pos[0] < current_pos[0] + len_motif:
            if pos[1] > current_pos[1]:
                current_pos = pos
                continue
        else:
            filtered_positions.append(
                '{}:{}'.format(current_pos[0], current_pos[1])
            )
            current_pos = pos
    return '\n'.join(filtered_positions)
