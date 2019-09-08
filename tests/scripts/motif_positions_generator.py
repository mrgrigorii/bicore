import argparse
import os
import random


SCRIPT_PATH = os.path.dirname(__file__)
OUTPUT_FILE_NAME = os.path.join(
    SCRIPT_PATH, os.pardir, 'motif_positions', 'test_motif_pos'
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-position', type=int)

    args = parser.parse_args()
    generate_motif_positions(args.max_position, random_true_threshold=0.95)


def generate_motif_positions(max_pos, random_true_threshold):
    lines = ['AHR_HUMAN.H11MO.0.B\n']
    
    for i in range(max_pos):
        random_value = random.random()
        if random_value > random_true_threshold:
            lines.append('{}:{}\n'.format(str(i), '1.0'))

    with open(OUTPUT_FILE_NAME, 'w') as out:
        out.write(''.join(lines))


if __name__ == '__main__':
    main()
