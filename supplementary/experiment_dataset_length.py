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

def get_nb_seqs_aver_len_group_families(dir_path):
    # 1) get all the files in the directory
    list_files = get_files_paths(dir_path)
    print("list files:", list_files)
    list_counts_len_avr = []
    for f in list_files:
        nb_seqs, avr_seq_len = count_and_average_len_seqs(f) 
        print("file {} = {} seqs , {} avr len".format(os.path.basename(f),nb_seqs,avr_seq_len))
        list_counts_len_avr.append((nb_seqs, avr_seq_len))

    return list_counts_len_avr

def scatter_list_tuples(list_nb_seqs_avr_len, test_name, is_log_scale = False):
    list_nb_seqs, list_avr_len = zip(*list_nb_seqs_avr_len)
    print(list_nb_seqs)
    print(list_avr_len)
    plt.scatter(list_nb_seqs, list_avr_len)
    plt.title(test_name)
    plt.xlabel("Number of sequences")
    plt.ylabel("Average sequences length")
    if is_log_scale:
        plt.xscale('log')
        plt.yscale('log')
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
    else:
        print (" dir_all_train_test_files exists already ...")
        # check if directory is empty
        nb_items = len(os.listdir(dir_all_train_test_files))
        if nb_items == 0:
            print (" dir_all_train_test_files is empty")
        else:
            print (" dir_all_train_test_files contains {} items".format(nb_items))  
        
        return

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

# note: noise data set, the file dosn't have exactly the same names
#       but the same prefix only.
#     TODO: write merge function for train test for that.
def merge_train_test_files_noise(main_directory, dir_out_name):
    dir_test = os.path.join(main_directory, "Test")
    dir_train = os.path.join(main_directory, "Train")
    dir_all_train_test_files = os.path.join(main_directory, dir_out_name)
    
    if not os.path.exists(dir_all_train_test_files):
        os.makedirs(dir_all_train_test_files)

    # since names are in RFxxxxx_nb-seqs example (RF00309_210),
    # RFxxxxx  are unique names
    # RFxxxxx  is the same in Train and Test files
    # we can get files of test and train in separate list, and sort theme
    # hence we will have the same files in the same order (because of uniqueness of RFxxxxx)
    # after that we can iterate by thier ids.
    # and creat a new file in all_train_test with RFxxxxx names (without the suffix of nb-seqs)

    list_files_paths_names_test = sorted(get_files_paths(dir_test))
    list_files_paths_names_train = sorted(get_files_paths(dir_train))
    
    for i in range(len(list_files_paths_names_test)):
        file_train = list_files_paths_names_train[i]
        file_test = list_files_paths_names_test[i]
        family_name = os.path.basename(file_train).split("_")[0] + ".fasta"# family name RFxxxxx
        file_train_test_out = os.path.join(dir_all_train_test_files, family_name)

        list_files_in = [file_train, file_test]
        with open(file_train_test_out, 'w') as outfile:
            for file_in in list_files_in:
                with open(file_in) as infile:
                    for line in infile:
                        outfile.write(line)


def main():
    dir_rfam_normal = r"C:\Users\ibra\Desktop\Infernal\nbF_all_nbSeqs_min_3"
    dir_noise = r"C:\Users\ibra\Desktop\Infernal\deep_ncrna_datasets\original"
    dir_clan_36 = r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_36_Train_Test"
    dir_clan_36_rfam = r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_36_Train_Test_Rfam"
    
    test_name_1 = "1st Dataset \'All Rfam\'"
    test_name_2 = "2nd Dataset \'Noise based\'"
    test_name_3 = "3rd Dataset \'36 Clans\'"
    test_name_4 = "4th Dataset \'Rfam family from 36 Clans\'"

    working_dir = dir_clan_36_rfam
    test_name = test_name_4

    dir_out_name = "all_train_test"
    #merge_train_test_files(working_dir, dir_out_name)
    merge_train_test_files_noise(working_dir, dir_out_name)
    dir_all_train_test_files = os.path.join(working_dir, dir_out_name)

    res = get_nb_seqs_aver_len_group_families(dir_all_train_test_files)

    print(" all results : ")
    print(res)
    #graph_list_tuples(res)
    scatter_list_tuples(res, test_name)

    # at the end : remove the all train test directory and its files.
    ##shutil.rmtree(dir_all_train_test_files)



if __name__ == '__main__':
    main()