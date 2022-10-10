"""
this file is only to use some script to faciltate some processing of my files
"""

# 1) Rename all files in the directory

import os
import sys
import shutil
import glob
import argparse
import subprocess
import time

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

def main():
    directory = (r"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_01-51-69_families\Train")
    #rename_all_files_in_directory(directory, "txt")
    rename_all_files_in_directory_prefix(directory, "red_")

if __name__ == "__main__":
    main()
