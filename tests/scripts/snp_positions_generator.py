import argparse
import os
import random


SCRIPT_PATH = os.path.dirname(__file__)
OUTPUT_FILE_NAME = os.path.join(
    SCRIPT_PATH, os.pardir, 'sorted_snp_positions_test'
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-position', type=int)

    args = parser.parse_args()
    generate_snp_positions(args.max_position, random_true_threshold=0.95)


def generate_snp_positions(max_pos, random_true_threshold):
    lines = []
    for i in range(max_pos):
        random_value = random.random()
        if random_value > random_true_threshold:
            lines.append('%s\n' % str(i))

    with open(OUTPUT_FILE_NAME, 'w') as out:
        out.write(''.join(lines))


if __name__ == '__main__':
    main()
