import gzip


def main():
    #snp_iter = snp_iteritems('simple_somatic_mutation.open(1).tsv.gz')
    snp_iter = snp_iteritems('C:/Games/Data/simple_somatic_mutation.open(1).tsv.gz')
    out = open('icgc_all_1chr_snp(2)', 'w')

    header = snp_iter.__next__()

    n = 0
    nn = 0
    data_buffer = ['\t'.join([header[1], header[2], header[8], header[9], header[17], header[18], header[19], header[20]])]
    for line in snp_iter:
        try:
            icgc_donor_id = line[1]
            project_code = line[2]
            chromosome = line[8]
            chromosome_start = line[9]
            quality_score = line[17]
            probability = line[17]
            total_read_count = line[17]
            mutant_allele_read_count = line[17]
        except:
            continue

        if chromosome == '1':
            n += 1
            data_buffer.append('\t'.join([
                icgc_donor_id,
                project_code,
                chromosome,
                chromosome_start,
                quality_score,
                probability,
                total_read_count,
                mutant_allele_read_count,
            ]))

            if n == 100000:
                print('Dumping', nn)
                nn += 1
                n = 0
                out.write('\n'.join(data_buffer))
                out.write('\n')
                data_buffer = []
    if n != 0:
        out.write('\n'.join(data_buffer))
        out.write('\n')


def snp_iteritems(file_path):
    with gzip.open(file_path, 'rb') as fin:
        data = fin.readline().decode('utf-8').split('\t')
        yield data
        while data:
            data = fin.readline().decode('utf-8').split('\t')
            yield data
        return None


if __name__ == '__main__':
    main()
