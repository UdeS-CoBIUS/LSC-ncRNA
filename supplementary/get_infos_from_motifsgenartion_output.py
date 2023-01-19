

""" the file is like : 
nb families : 10
Min nb seqs : 20
Max nb seqs : 40
Min length motif : 12
Max length motif : 13
 False 0
BETA : 0
Alpha : -1
nbOccrs_allowed : 1
 Total nb_seqs = 6366
 Total nb_motifs = 78818
 sum_cm_nb_occrs_alpha = 0
 Average nb_occres = 0
 Total nb_seqs = 6366
 Total nb_motifs = 78818
file_csv.close();
Time taken before save_to_csv is : 7.220197 sec
Time taken by whole program is : 45.896144983 sec
"""
# we want from it : 
# Total nb_seqs = 6366
# Total nb_motifs = 78818
# Time taken before save_to_csv is : 7.220197 sec
# Time taken by whole program is : 45.896144983 sec
# and return it as :
# nb_seqs, nb_motifs, time_generated, all time
# 6366, 78818 , 7.220197 , 45.896144983


import os
import csv
import sys

def get_infos_results_output(filename_path):

    nb_seqs, nb_motifs, time_generated, time_generation_save_csv = 0, 0, 0, 0

    with open(filename_path, "r") as in_file:
        
        for line in in_file:
            # check for nb_seqs
            if " Total nb_seqs" in line:
                nb_seqs = line.split("=")[1].strip()
            # check for nb_motifs
            if " Total nb_motifs" in line:
                nb_motifs = line.split("=")[1].strip()
            # check for time_generated
            if "Time taken before save_to_csv is" in line:
                time_generated = line.split("is :")[1].strip().split()[0]
            # check for all_time
            if "Time taken by whole program is" in line:
                time_generation_save_csv = line.split("is :")[1].strip().split()[0]
    
    file_name_basename = os.path.basename(filename_path)
    file_name_only = os.path.splitext(file_name_basename)[0]

    return file_name_only, nb_seqs, nb_motifs, time_generated, time_generation_save_csv

def get_files_paths(directory, prefix="out_F_"):
    files = [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and os.path.basename(f).startswith(prefix)]
    return files


def get_infos_dir(directory):
    
    list_files_paths = get_files_paths(directory)

    all_results = []

    for file in list_files_paths:
        result = get_infos_results_output(file)
        all_results.append(result)
        print(result)

    # sort results according to file names
    sorted(all_results, key=lambda x: x[0])
    # write results to csv file
    with open("all_results_out_F_.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(all_results)
    
    print("All results out written to csv file [out_F_all_results.csv]")


def main():
    my_dir = sys.argv[1]
    print(" dir = " + my_dir)
    get_infos_dir(my_dir)


if __name__ == '__main__':
    main()
