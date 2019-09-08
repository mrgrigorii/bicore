import gzip
import os

lines_in_vcf = 6468094

def load_fasta_from_gzip(fasta_gzip_name):
    with gzip.open(fasta_gzip_name, 'rb') as f:
        data = f.read()
        return load_fasta_data(data)
        

def load_fasta_data(data):
    result_data = []
    data = data.strip('>').split('>')
    print('Fasta file have ' + str(len(data)) + ' scaffolds')
    for scaffold in data:
        scaffold = scaffold.strip('\n').split('\n')
        name = scaffold[0]
        scaffold = ''.join(scaffold[1:])
        result_data.append([name, scaffold])
    data = None
    return result_data


def load_hocomoco_motif_data(file_name):
    with open(file_name, 'r') as f:
        data = f.read().strip('\n').split('\n')
        motif_name = data[0]
        data = [i.split(':') for i in data[1:]]
        data = [[int(i[0]), float(i[1])] for i in data]
    return motif_name, data

def get_motif_len(motif_name):
    with open('HOCOMOCOv11_core_pwms_HUMAN_mono.txt', 'r') as f:
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

def load_gen_pos(file_name):
    list_pos = []
    num = 0
    with open(file_name, 'r') as f:
        for line in f:
            if line[0:4] == 'chr1':
                line = line.split()
                if line[6] == '+':
                    list_pos.append([int(line[3]), int(line[4])])
                    num += 1
            if num >= 2000:
                return list_pos
    return list_pos
                

#gen_pos_list = load_gen_pos('gencode.v28.basic.annotation.gtf')
#gen_pos_list.sort()
#print(len(gen_pos_list))
#input()

#chrom_1 = load_fasta_from_gzip('asp_2/chr1.fna.gz')
#print(len(chrom_1[0][1]))
#input()

motif_data_list = []
result_list = []
motif_files = os.listdir('dataD/')
print(motif_files)

for motif in motif_files:
    motif_data = load_hocomoco_motif_data('dataD/' + motif)
    list_pos = get_list_pos_motif(motif_data)
    motif_data_list.append([motif, list_pos])
    
    header = '#' + motif + '|||num_positions:' + str(len(list_pos))
    result_list.append([header])

#GENE
#list_pos = gen_pos_list
out = open('test.vcf', 'w')

file_name = 'chr_1.txt.gz'
n = 0
with gzip.open(file_name, 'r') as f:
    N = 0
    while True:
        lines = f.readlines(1000)
        if not lines:
            break

        N2 = 0
        start = False
        for line in lines:
            line = line.decode()
            print(line, N, N2)
            N += 1
            N2 += 1
            if N >= 5000:
                out.close()
                input('Stop')
            out.write(line)
            line = line.split('\t')
            if len(line) > 10 and line[6] == '1': # 1 chromosome hints
                n += 1
                if n % 10000 == 0:
                    print(str(n*100.0/lines_in_vcf) + '%')
                pos = int(line[11])
                print('pos', pos, line)
                #print(len(chrom_1), len(chrom_1[0][1]), pos)
                #print(chrom, pos, ref, alt, info)
                #print(chrom, pos, ref, alt, chrom_1[0][1][pos-1])
                
                
                pos_2 = pos #- 1   # my_pos = real_pos - 1
                
                
                for i_in_res, motif in enumerate(motif_data_list):

                    remove_list = []
                    for i in motif[1]:
                        if i[0]<=pos_2<=i[1]:
                            chrom = line[0]
                            snp_id = line[2]
                            ref = line[3]
                            alt = line[4]
                            info = line[7]
                            result_list[i_in_res].append('\t'.join(map(str,[chrom, pos, ref, alt, info, i])))
                            #print('DAAA', pos)
                        elif i[1] < pos:
                            if i[1] + 1000 < pos_2:
                                remove_list.append(i)
                            continue
                        else: 
                            break
                    
                    for i in remove_list:
                        motif[1].remove(i)
                
                #if len(list_pos) == 0:
                    #print('Break')
                    #break

    print('Lines: ' + str(n))
                    
        #n += 1
        #if n == 10000:
            #break

for i, motif in enumerate(result_list):
    out = open(motif_data_list[i][0] + '.txt', 'w')
    out.write('\n'.join(motif))
    out.close()
