import os

"""
this code is for merging files with same names that are in two different directories into a file with same name into output directory
this code is with help of chatGPT generator
"""

# set directory paths
input_dir1 = "/path/to/dir1"
input_dir2 = "/path/to/dir2"
output_dir3 = "/path/to/dir3"

# get list of files in dir1 and dir2
files1 = os.listdir(input_dir1)
files2 = os.listdir(input_dir2)

# Iterate through the list of file names and merge the files with the same name
file_nb = 0
for file_name in input_dir1:
    file_nb = +1
    print(" {} : {}".format(file_nb, file_name))
    if file_name in input_dir2:  # if file with the same name is in dir2
        file1_path = os.path.join(input_dir1, file_name)
        file2_path = os.path.join(input_dir2, file_name)
        file3_path = os.path.join(output_dir3, file_name)
        with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2, open(file3_path, 'w') as file3:
            file3.write(file1.read())
            file3.write(file2.read())

print "finished :) "