import os
import time


def test_runtime(func):
    start_time = time.time()

    def wrupped_func(*args, **kwargs):
        func(*args, **kwargs)
        runtime = round(time.time() - start_time, 2)
        print('Runtime: %s sec' % (str(runtime)))

    return wrupped_func


def load_hocomoco_motif(file_name):
    motif_list = []
    with open(file_name, 'r') as f:
        data = f.read().split('>')[1:]
        for motif in data:
            motif = motif.split('\n')
            motif_name = motif[0]
            motif_matrix = list(
                list(map(float, i.split())) + [0.0] for i in motif[1:] if i
            )
            motif_list.append([motif_name, motif_matrix])

    return motif_list


def best_motif_score(motif):
    best_score = 0
    for i in motif[1]:
        best_score += max(i)
    return best_score


def worst_motif_score(motif):
    worst_score = 0
    for i in motif[1]:
        worst_score += min(i)
    return worst_score


def dump_data(data, filename):
    check_folder_exists(filename)
    with open(filename, 'w') as out:
        out.write(data)


def check_folder_exists(filename):
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
