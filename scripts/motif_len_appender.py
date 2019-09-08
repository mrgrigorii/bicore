import os


def main():
    motif_files = os.listdir('data/')

    for motif in motif_files:
        with open('data/' + motif, 'r') as fin:
            data = fin.read().strip('\n').split('\n')
            motif_name = data[0]
            len_motif = get_motif_len(motif_name)
            data[0] += ':' + str(len_motif)

            with open('updated_results/' + motif, 'w') as out:
                out.write('\n'.join(data))


def get_motif_len(motif_name):
    with open('HOCOMOCOv11_core_pwms_HUMAN_mono.txt', 'r') as f:
        data = f.read().split('>')[1:]
        for motif in data:
            if motif_name in motif:
                return len(motif.strip('\n').split('\n')) - 1


if __name__ == "__main__":
    main()
