"""
this file is only to use some script to faciltate some processing of my files
"""

# 1) Rename all files in the directory

import requests
from bs4 import BeautifulSoup


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
#<span class="listItem">
#      <a title="RF00017" href="https://rfam.xfam.org/family/RF00017">
#        Metazoa_SRP</a>
#    </span>
# we need to get items from class="listItem", get the title "RF00017"
# we need to get all items from class="listItem".
def get_clan_Rfam_classes():

    list_rfam = []

    r = requests.get("https://rfam.xfam.org/clan/CL00003")
    soup = BeautifulSoup(r.text, "html.parser")
    list_item = soup.find_all(class_="listItem")
    for item in list_item:
        list_rfam.append(item.find("a").get("title"))

    return list_rfam
    
    

def main():
    #directory = (r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_01-51-69_families\Train")
    #rename_all_files_in_directory(directory, "txt")
    #rename_all_files_in_directory_prefix(directory, "red_")
    print("hello word families")
    list_rfam = get_clan_Rfam_classes()
    print(list_rfam)

if __name__ == "__main__":
    main()
