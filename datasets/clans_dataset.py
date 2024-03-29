"""
this file is only to use some script to faciltate some processing of my files
"""

# 1) Rename all files in the directory

import requests
from bs4 import BeautifulSoup
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import shutil

max_nb_seqs = 12  # max numbers of sequence by fasta file.


# get all files path that begin with some prefix
def get_all_files_path(directory, prefix):
    all_files_path = []
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.startswith(prefix):
                all_files_path.append(os.path.join(root, f))
    return all_files_path


# get all files path that end with ".txt"
def get_files_path(directory, extension):
    files_path = []
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.endswith(extension):
                files_path.append(os.path.join(root, f))
    return files_path


# 1) Rename all files in the directory

def rename_all_files_in_directory(directory, extension):
    file_path = get_files_path(directory, extension)

    for f in file_path:
        if os.path.isfile(f):
            # get the file name without all extensions
            # get the base name of the file
            # split the name and get the first name only.
            base_name = os.path.splitext(os.path.basename(f))[0]
            file_name = base_name.split('.')[0]
            new_file_name = os.path.join(directory, "{}.fasta".format(file_name))
            # rename the file to the new_file_name
            os.rename(f, new_file_name)


def rename_all_files_in_directory_prefix(directory, prefix="red_"):
    file_path = get_all_files_path(directory, prefix)

    for f in file_path:
        if os.path.isfile(f):
            # split the name using "_"
            base_name = os.path.splitext(os.path.basename(f))[0]
            file_name = base_name.split('_')[1]
            new_file_name = os.path.join(directory, "{}.fasta".format(file_name))
            # rename the file to the new_file_name
            os.rename(f, new_file_name)


# get all RFam classes from a given Clan, by its html link.
# urls: https://rfam.xfam.org/clan/CL00003
# we have an html tag like this:
# <span class="listItem">
#      <a title="RF00017" href="https://rfam.xfam.org/family/RF00017">
#        Metazoa_SRP</a>
#    </span>
# we need to get items from class="listItem", get the title "RF00017"
# we need to get all items from class="listItem".
def get_clan_Rfam_classes_names(clan_url):
    list_rfam = []

    r = requests.get(clan_url)
    soup = BeautifulSoup(r.text, "html.parser")
    list_item = soup.find_all(class_="listItem")
    for item in list_item:
        list_rfam.append(item.find("a").get("title"))

    return list_rfam


def get_clans_names(url_clans_page):
    list_clans_names = []

    r = requests.get(url_clans_page)
    soup = BeautifulSoup(r.text, "html.parser")
    # list_item = soup.find_all(class_="listItem")
    # for item in list_item:
    #    list_rfam.append(item.find("a").get("title"))
    # todo: to be completed later.

    return list_clans_names


# download file from urls.
def download_file(url, file_path_name):
    r = requests.get(url)
    with open(file_path_name, 'wb') as f:
        f.write(r.content)


def get_download_clan_Rfam_classes(url_clan, directory):
    clan_file_name = url_clan.split("/")[-1]
    print(" clan_file_name: " + clan_file_name)

    # create a folder for "clan_file_name" in the directory
    dir_clan = os.path.join(directory, clan_file_name)

    # check if folder exists, stop the program and return.
    if os.path.exists(dir_clan):
        print("The folder {} does exist... Exiting...".format(dir_clan))
        return

    print(" creating folder: " + dir_clan)
    os.makedirs(dir_clan)
    print(" creating folder: " + dir_clan)

    # get all Rfam classes names from the url_clan
    list_rfam = get_clan_Rfam_classes_names(url_clan)
    print(" list_rfam: " + str(list_rfam))

    # link of RFam downlowd alignments sequen, ungraped fasta file
    # https://rfam.xfam.org/family/RF00169/alignment?acc=RF00169&format=fastau&download=1

    # download all RFam classes
    for rfam in list_rfam:
        file_path_name = os.path.join(dir_clan, rfam + ".fasta")
        url_rfam_align_ungped = "https://rfam.xfam.org/family/{}/alignment?acc={}&format=fastau&download=1".format(rfam,
                                                                                                                   rfam)
        print("file_path_name={}".format(file_path_name))
        print("url_rfam_align_ungped={}".format(url_rfam_align_ungped))
        download_file(url_rfam_align_ungped, file_path_name)
        print(" finished download_file----------------------")


def clans_family_download(directory, base_url_clan, list_clan_names):
    for clan_name in list_clan_names:
        url_clan = base_url_clan + clan_name
        get_download_clan_Rfam_classes(url_clan, directory)


def sequences_count(fasta_file):
    nb_seqs = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                nb_seqs += 1
    return nb_seqs


# for each family in the folder, count the number of fasta sequences
# if there are more than 12 sequences, save only the 12 sequences, and remove the others.
def reduce_nb_sequences(fasta_file_in, fasta_file_out, max_nb_sequences):
    print("fasta_file_in: ", fasta_file_in)
    print("fasta_file_out: ", fasta_file_out)

    idx = 0
    with open(fasta_file_in, 'r') as fasta_in:
        with open(fasta_file_out, 'w') as fasta_out:
            for line in fasta_in:
                if line.startswith('>'):
                    print(" seq : {}".format(idx))
                    if idx <= max_nb_sequences:
                        print(" write seq in file ")
                        fasta_out.write(line)
                        idx += 1
                    else:
                        break  # finish the writing.
                else:
                    fasta_out.write(line)
    print(" finished writing")


def reduce_nb_sequences_folder(folder, nb_sequences):
    # get all fasta file in the directory
    fasta_files = get_files_path(folder, "fasta")
    print(" nb fasta files: ", len(fasta_files))

    for fasta_file in fasta_files:
        fasta_file_name = os.path.basename(fasta_file)
        red_fasta_file_name = os.path.join(os.path.dirname(fasta_file), "red_" + fasta_file_name)
        print(" red_fasta_file_name: ", red_fasta_file_name)
        reduce_nb_sequences(fasta_file, red_fasta_file_name, nb_sequences)


def remove_files_bigin_prefix(folder, prefix):
    # get all fasta file in the directory
    fasta_files = get_files_path(folder, "fasta")
    for filename in fasta_files:
        fasta_file_name = os.path.basename(filename)
        if fasta_file_name.startswith(prefix):
            os.remove(filename)


def rename_files_remove_prefix_from_name(folder, prefix):
    # get all fasta file in the directory
    fasta_files = get_files_path(folder, "fasta")
    for filename in fasta_files:
        fasta_file_name = os.path.basename(filename)
        new_name = fasta_file_name[len(prefix):]
        if fasta_file_name.startswith(prefix):
            path_new_name = os.path.join(folder, new_name)
            os.rename(filename, path_new_name)


def clans_family_reducer(directory, list_clan_names):
    print(" In clans_family_reducer ................")

    for clan_name in list_clan_names:
        url_clan = directory + clan_name
        print(" Reduce clan : " + clan_name)
        reduce_nb_sequences_folder(url_clan, max_nb_seqs)
        remove_files_bigin_prefix(url_clan, "RF")
        rename_files_remove_prefix_from_name(url_clan, "red_")


def generate_train_test_files(folder):
    # 1) create directory for train test
    folder_name = os.path.basename(folder)
    folder_dir = os.path.dirname(folder)
    train_test_folder = os.path.join(folder_dir, folder_name + "_train_test")
    if not os.path.exists(train_test_folder):
        os.makedirs(train_test_folder)

    # 2) lunch
    percentage_nb_seqs_train = 70
    construct_train_test_files(folder, train_test_folder, percentage_nb_seqs_train)


def construct_train_test_files(path_dir_in, path_dir_out, percentage_nb_seqs_train):
    # 1) create Train and Test folder
    train_folder = os.path.join(path_dir_out, "Train")
    test_folder = os.path.join(path_dir_out, "Test")
    res_folder = os.path.join(path_dir_out, "res")  # we need this folder during experiment

    if not os.path.exists(train_folder):
        os.makedirs(train_folder)
    if not os.path.exists(test_folder):
        os.makedirs(test_folder)
    if not os.path.exists(res_folder):
        os.makedirs(res_folder)

    # get all fasta files in path_dir
    fasta_files = get_files_path(path_dir_in, "fasta")
    print(" nb fasta files: ", len(fasta_files))

    for fasta_file in fasta_files:
        fasta_file_name = os.path.basename(fasta_file)
        file_train = os.path.join(train_folder, fasta_file_name)
        file_test = os.path.join(test_folder, fasta_file_name)
        nb_seqs_file = sequences_count(fasta_file)
        nb_seqs_train = (nb_seqs_file * percentage_nb_seqs_train) / 100
        idx = 0
        train_records = []
        test_records = []

        for record in SeqIO.parse(fasta_file, "fasta"):
            if idx < nb_seqs_train:
                train_records.append(record)
                idx += 1
            else:
                test_records.append(record)

        SeqIO.write(train_records, file_train, "fasta")
        SeqIO.write(test_records, file_test, "fasta")


def clans_train_test(directory, list_clan_names):
    for clan_name in list_clan_names:
        url_clan = directory + clan_name
        generate_train_test_files(url_clan)


def gather_sequences_multiples_files_into_one_file(dir_clan_name, dir_files_in, path_file_out):
    # get all fasta files in path_dir
    list_fasta_files = get_files_path(dir_files_in, "fasta")
    print(" nb fasta files: ", len(list_fasta_files))

    clan_records = []

    idx = 0

    for fasta_file in list_fasta_files:
        rfam_file_name = os.path.basename(fasta_file)
        for record in SeqIO.parse(fasta_file, "fasta"):
            rec = SeqRecord(record.seq, id="{}_{}".format(dir_clan_name, idx),
                            description=" | {}".format(rfam_file_name))
            clan_records.append(rec)
            idx += 1

    # at the end , write the clan records of all families.
    SeqIO.write(clan_records, path_file_out, "fasta")


def gather_sequences_clan_train_test(clan_name, clan_dir, dir_out):
    # 1) create Train and Test folder
    dir_out_train = os.path.join(dir_out, "Train")  # train and test folder supposed created before.
    dir_out_test = os.path.join(dir_out, "Test")

    train_out_file_path_name = os.path.join(dir_out_train, "{}.fasta".format(clan_name))
    test_out_file_path_name = os.path.join(dir_out_test, "{}.fasta".format(clan_name))

    dir_in_train = os.path.join(clan_dir, "Train")
    dir_in_test = os.path.join(clan_dir, "Test")

    gather_sequences_multiples_files_into_one_file(clan_name, dir_in_train, train_out_file_path_name)
    gather_sequences_multiples_files_into_one_file(clan_name, dir_in_test, test_out_file_path_name)


def gather_train_test_families_in_Clans(directory, list_clan_names, dir_clan_out):
    for clan_name in list_clan_names:
        print(clan_name + " -------------- ")
        dir_clan_in = directory + "{}_train_test".format(clan_name)
        gather_sequences_clan_train_test(clan_name, dir_clan_in, dir_clan_out)


def clan(main_directory, dir_all_clans_out, base_url_clan, list_clan_names):
    """
    :param main_directory: main directory where we downlad clan temporray families, and work on them.
    :param dir_all_clans_out: folder where we gather all sequnces of clans in Train and Test sub folder.
    :param base_url_clan:  the clan url (example: https://rfam.xfam.org/clan/CL00051)
    :param list_clan_names: list names of clans that will be treated.
    :return: None
    1) download all clan family is separate temporary folder
    2) reduce the number of sequences  in each fasta file by max_nb_seqs
    3) split each file is Train Test part.
    4) gather all train test files of individual rfam family in their main clan faile Train Test.
    """
    clans_family_download(main_directory, base_url_clan, list_clan_names)
    clans_family_reducer(main_directory, list_clan_names)
    clans_train_test(main_directory, list_clan_names)

    gather_train_test_families_in_Clans(main_directory, list_clan_names, dir_all_clans_out)


def copy_files(source_dir, target_dir):
    # get all files fro source dir:
    list_files = get_files_path(source_dir, "fasta")
    for file in list_files:
        shutil.copy(file, target_dir)


def gather_train_test_rfam(main_directory, rfam_only_dir, list_clan_names):
    test_dir_out = os.path.join(rfam_only_dir, "Test") # the sub folder supposed to be created manually
    train_dir_out = os.path.join(rfam_only_dir, "Train")

    for clan_name in list_clan_names:
        clan_dir = os.path.join(main_directory, clan_name + "_train_test")
        clan_dir_train = os.path.join(clan_dir, "Train")
        clan_dir_test = os.path.join(clan_dir, "Test")

        copy_files(clan_dir_train, train_dir_out)
        copy_files(clan_dir_test, test_dir_out)


# replace by creating a second output file, then remove the first one, and rename the second one.
def replace_seqs_id_by_family_id_file(file_path):
    file_base_name = os.path.basename(file_path)
    file_dir_path = os.path.dirname(file_path)
    file_out_base_name = file_base_name + "_out"
    file_out_dir_path = os.path.join(file_dir_path, file_out_base_name)

    # 1) write all line to out, with changing the line id start with ">"
    idx = 0
    with open(file_path, 'r') as fasta_in:
        with open(file_out_dir_path, 'w') as fasta_out:
            for line in fasta_in:
                if line.startswith('>'):
                    fasta_out.write(">{}_{}\n".format(file_base_name, idx))
                    idx += 1
                else:
                    fasta_out.write(line)
    # 2) remove original file
    os.remove(file_path)

    # 3) rename out file to original.
    os.rename(file_out_dir_path, file_path)


def replace_seqs_id_by_family_id_dir(directory):
    # 1) read all files in directory
    list_files = get_files_path(directory, "fasta")

    for file in list_files:
        replace_seqs_id_by_family_id_file(file)


def rfam_individual_family_only_from_clans(main_directory, rfam_only_dir, list_clan_names):
    """
    :param main_directory: main directory where we gathered all clans and their rfam families (Train and Test)
    :param rfam_only_dir: folder that contain (Train, Test) rfam only files.
    :param list_clan_names: list of used clans names
    :return: None

    1) copy individual rfam train and test files from the clan folder into one folder "rfam_only_dir" (train and test)
    2) replace seqs id by family id in rfam_only_dir/Train and Test rfam files.
    """
    # 1:
    gather_train_test_rfam(main_directory, rfam_only_dir, list_clan_names)

    # 2:
    dir_rfam_only_train = os.path.join(rfam_only_dir, "Train")
    dir_rfam_only_test = os.path.join(rfam_only_dir, "Test")
    replace_seqs_id_by_family_id_dir(dir_rfam_only_train)
    replace_seqs_id_by_family_id_dir(dir_rfam_only_test)


def main():

    mode = "rf" # cl: clans,  rf: rfam form clans, clrf: cl+rf

    directory = r"C:/Users/ibra/Desktop/Infernal/Clans ncRNA/"
    dir_all_clan_out = r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_36_Train_Test"  # created manualy, and also Train, Test, res sub folder
    rfam_only_dir = r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_36_Train_Test_Rfam"
    base_url_clan = "https://rfam.xfam.org/clan/"
    list_clan_names = ["CL00051",
                       "CL00003",
                       "CL00069",
                       "CL00106",
                       "CL00038",
                       "CL00054",
                       "CL00002",
                       "CL00014",
                       "CL00102",
                       "CL00001",
                       "CL00117",
                       "CL00021",
                       "CL00057",
                       "CL00111",
                       "CL00118",
                       "CL00066",
                       "CL00112",
                       "CL00116",
                       "CL00010",
                       "CL00005",
                       "CL00027",
                       "CL00053",
                       "CL00063",
                       "CL00100",
                       "CL00119",
                       "CL00004",
                       "CL00032",
                       "CL00034",
                       "CL00035",
                       "CL00040",
                       "CL00096",
                       "CL00093",
                       "CL00015",
                       "CL00121",
                       "CL00036",
                       "CL00045"
                       ]

    if mode == "cl":
        clan(directory, dir_all_clan_out, base_url_clan, list_clan_names)
    elif mode == "rf":
        rfam_individual_family_only_from_clans(directory, rfam_only_dir, list_clan_names)
    else: # both, or we can check "clrf"
        clan(directory, dir_all_clan_out, base_url_clan, list_clan_names)
        rfam_individual_family_only_from_clans(directory, rfam_only_dir, list_clan_names)





if __name__ == "__main__":
    main()
