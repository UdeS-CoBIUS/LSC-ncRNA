# experiment_dataset_length:
#  generate a graph where on the X axis put family size (number of seqs),
#  and on Y axis the average len of seqs.
#  this for 3 dataset:
#  1) Rfam dataset,
#  2) Noise dataset for 0 noise,
#  3) Clan dataset (and for Rfam from clans)
#
# 	  for that we have compute for each family 1) number of seqs, and the average len seqs.
import os
import matplotlib.pyplot as plt
from Bio.SeqIO.FastaIO import SimpleFastaParser
import shutil

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
    
    count = 0
    total_len = 0
    with open(file_family) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            count += 1
            total_len += len(seq)
            #print(" seq num =", count, " | len =",len(seq))
            
    return count, total_len//count # integer division

def get_files_paths(directory):
    files = [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    return files

def group_families(dir_path):
    # 1) get all the files in the directory
    list_files = get_files_paths(dir_path)
    print("list files:", list_files)
    list_counts_len_avr = []
    for f in list_files:
        nb_seqs, avr_seq_len = count_and_average_len_seqs(f) 
        print("file {} = {} seqs , {} avr len".format(os.path.basename(f),nb_seqs,avr_seq_len))
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


def merge_train_test_files(main_directory, dir_out_name):
    dir_test = os.path.join(main_directory, "Test")
    dir_train = os.path.join(main_directory, "Train")
    dir_all_train_test_files = os.path.join(main_directory, dir_out_name)
    
    if not os.path.exists(dir_all_train_test_files):
        os.makedirs(dir_all_train_test_files)

    list_files_names = [f for f in os.listdir(dir_train) if os.path.isfile(os.path.join(dir_train, f))] # the files have the same names in Test and Train

    
    for file_name in list_files_names:
        file_train = os.path.join(dir_train, file_name)
        file_test = os.path.join(dir_test, file_name)
        file_train_test_out = os.path.join(dir_all_train_test_files, file_name)

        list_files_in = [file_train, file_test]
        with open(file_train_test_out, 'w') as outfile:
            for file_in in list_files_in:
                with open(file_in) as infile:
                    for line in infile:
                        outfile.write(line)


def main():
    
    dir_clan_36 = r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_36_Train_Test"
    dir_clan_36_rfam = r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_36_Train_Test_Rfam"

    dir_out_name = "all_train_test"
    merge_train_test_files(dir_clan_36_rfam, dir_out_name)
    dir_all_train_test_files = os.path.join(dir_clan_36_rfam, dir_out_name)

    res = group_families(dir_all_train_test_files)

    print(" all results : ")
    print(res)
    #graph_list_tuples(res)
    scatter_list_tuples(res)

    # at the end : remove the all train test directory and its files.
    shutil.rmtree(dir_all_train_test_files)



if __name__ == '__main__':
    main()