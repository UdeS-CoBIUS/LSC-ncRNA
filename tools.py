"""
this file is only to use some script to faciltate some processing of my files
"""

# 1) Rename all files in the directory

import requests
from bs4 import BeautifulSoup
import os


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
        print("The folder does exist.")
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
        url_rfam_align_ungped = "https://rfam.xfam.org/family/{}/alignment?acc={}&format=fastau&download=1".format(rfam, rfam)
        print("file_path_name={}".format(file_path_name))
        print("url_rfam_align_ungped={}".format(url_rfam_align_ungped))
        download_file(url_rfam_align_ungped, file_path_name)
        print(" finished download_file----------------------")



def main():
    # directory = (r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_01-51-69_families\Train")
    # rename_all_files_in_directory(directory, "txt")
    # rename_all_files_in_directory_prefix(directory, "red_")

    directory = r"C:/Users/ibra/Desktop/Infernal/Clans ncRNA/"
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
                       "CL00111",
                       "CL00118",
                       "CL00066",
                       "CL00112",
                       "CL00116",
                       "CL00010",
                       "CL00005",
                       "CL00027"
                       ]
    for clan_name in list_clan_names:
        url_clan = base_url_clan + clan_name
        get_download_clan_Rfam_classes(url_clan, directory)


if __name__ == "__main__":
    main()
