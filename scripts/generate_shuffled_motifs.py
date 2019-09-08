import random

from bicore import utils


def main():
    random_motifs = _get_shuffled_motifs()
    _save_motifs(random_motifs, 'shuffled_motifs_pwms')


def _get_shuffled_motifs():
    motifs = utils.load_hocomoco_motif('HOCOMOCOv11_core_pwms_HUMAN_mono.txt')

    suffled_motifs = []
    for motif in motifs:
        name = motif[0]
        nodes = motif[1]
        random.shuffle(nodes)
        new_name = name + '_shuffled'
        suffled_motifs.append([new_name, nodes])
    return suffled_motifs


def _save_motifs(motifs, filename):
    with open(filename, 'w') as out:
        for motif in motifs:
            out.write('>%s\n' % motif[0])
            for node in motif[1]:
                out.write('\t'.join(list(map(str, node))) + '\n')


if __name__ == '__main__':
    main()
