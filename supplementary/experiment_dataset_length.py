# experiment_dataset_length:
#  generate a graph where on the X axis put family size (number of seqs),
#  and on Y axis the average len of seqs.
#  this for 3 dataset:
#  1) Rfam dataset,
#  2) Noise dataset for 0 noise,
#  3) Clan dataset (and for Rfam from clans)
#
# 	  for that we have compute for each family 1) number of seqs, and the average len seqs.

from os import listdir
from os.path import isfile, join, basename
import matplotlib.pyplot as plt

# count the number of sequences in each family
# with simply counting the sign of sequence id ">"
# we can use : `grep -c ">" file.fasta` but in shell, with `subprocess`
def count_number_of_seqs(file_family):
    nb_seqs = 0
    with open(file_family) as f:
        for line in f:
            if line.startswith('>'):
                nb_seqs += 1

    return nb_seqs

# this form : http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:low-level-fasta-fastq
def count_and_average_len_seqs(file_family):
    
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    count = 0
    total_len = 0
    with open(file_family) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            count += 1
            total_len += len(seq)
            #print(" seq num =", count, " | len =",len(seq))
            
    return count, total_len//count # integer division

def get_files_paths(directory):
    files = [join(directory, f) for f in listdir(directory) if isfile(join(directory, f))]
    return files

def group_families(dir_path):
    # 1) get all the files in the directory
    list_files = get_files_paths(dir_path)
    print("list files:", list_files)
    list_counts_len_avr = []
    for f in list_files:
        nb_seqs, avr_seq_len = count_and_average_len_seqs(f) 
        print("file {} = {} seqs , {} avr len".format(basename(f),nb_seqs,avr_seq_len))
        list_counts_len_avr.append((nb_seqs, avr_seq_len))

    return list_counts_len_avr

def scatter_list_tuples(list_nb_seqs_avr_len):
    list_nb_seqs, list_avr_len = zip(*list_nb_seqs_avr_len)
    print(list_nb_seqs)
    print(list_avr_len)
    plt.scatter(list_nb_seqs, list_avr_len)
    plt.show()

def graph_list_tuples(list_tuples):

    x_val = [x[0] for x in list_tuples]
    y_val = [x[1] for x in list_tuples]

    print(x_val)
    plt.plot(x_val,y_val)
    plt.plot(x_val,y_val,'or')
    plt.show()

def main():
    dir_test = r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\CL00001"
    dir_noise_0 = r"C:\Users\ibra\Desktop\Infernal\deep_ncrna_datasets\original\Train"
    dir_all_rfam_min_3_train_only = r"C:\Users\ibra\Desktop\Infernal\nbF_all_nbSeqs_min_3\Train"
    dir_clan_36 = r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_36_Train_Test_Rfam\Train"

    dir_test = dir_all_rfam_min_3_train_only
    res = group_families(dir_test)

    print(" all results : ")
    print(res)
    #graph_list_tuples(res)
    scatter_list_tuples(res)



if __name__ == '__main__':
    main()