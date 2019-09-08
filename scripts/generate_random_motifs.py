import random

from bicore import utils


def main():
    maxmin = _get_param_stats()
    random_motifs = _generete_test_motifs(num=200, maxmin=maxmin)
    _save_motifs(random_motifs, 'test_motifs_pwms_200')


def _get_param_stats():
    motifs = utils.load_hocomoco_motif('HOCOMOCOv11_core_pwms_HUMAN_mono.txt')
    max_score = 0
    min_score = 0
    for motif in motifs:
        for node in motif[1]:
            max_s = max(node)
            if max_s > max_score:
                max_score = max_s

            min_s = min(node)
            if min_s < min_score:
                min_score = min_s

    return min_score, max_score


def _generete_test_motifs(num=10, motif_len=[7, 15], maxmin=[0, 1]):
    random_motifs = []
    for i in range(num):
        random_motif_name = 'random_motif_%s' % str(i)
        random_motif_len = random.randint(*motif_len)

        random_motif = []
        for j in range(random_motif_len):
            line = _generate_line(maxmin)
            random_motif.append(line)

        random_motifs.append([random_motif_name, random_motif])

    return random_motifs


def _generate_line(maxmin):
    good_line = False
    while not good_line:
        good_line = False
        line = [random.uniform(*maxmin) for _ in range(4)]
        for n in line:
            if n > 0.0:
                good_line = True
    return line


def _save_motifs(motifs, filename):
    with open(filename, 'w') as out:
        for motif in motifs:
            out.write('>%s\n' % motif[0])
            for node in motif[1]:
                out.write('\t'.join(list(map(str, node))) + '\n')


if __name__ == '__main__':
    main()
