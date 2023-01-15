
import math
import csv

# code help chatGPT :)

"""
given an initial list of numbers L=[nb1, nb2,...,nbn], with a total of sum(L)
we need to generate a new list with a new total number N,
and we need to keep the same percenatge of the initial numbers.

The list must contain only integers
for that we use floor, and after that we begin to add 1 to the values inside the list till we complete the total N
in our case we begin to add 1 from the last element, because our initial list is ordered from greatest to least, and like this it may that the last value be 0 when scaling (reducing)

@N : the total number we need to scale for
@original_list: the original list
@returns the new list with a new total number N, and the numbers in it, respect (approximatly) the initial percenatge
"""
def scale_list_int(N, original_list):
    total = sum(original_list)
    ratio = N/total
    new_list = [math.floor(value * ratio) for value in original_list]
    diff = N - sum(new_list)
    if diff != 0:
        for i in range(diff):
            new_list[-1-i] += 1
    return new_list


def write_to_csv_dict_key_list_values(my_dict, csv_file_name):
    
    with open(csv_file_name, mode='w', newline='') as csv_file:
        fieldnames = list(my_dict.keys())
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()
        rows = [dict(zip(my_dict.keys(), row)) for row in zip(*my_dict.values())]
        for row in rows:
            writer.writerow(row)
        
        print(" finished writing: ",csv_file_name)



original_list = [1064, 500, 407, 191, 119, 4, 4]
list_totals_dataset_600 = [600,500,400,300,200,100]
list_totals_dataset_350 = [350,250,150,50]
list_new_totals_dataset = list_totals_dataset_350 + list_totals_dataset_600


def main():
    
    csv_file_name = "dataset_distribution_nb_seqs.csv"
    my_dict = {}

    for new_total in list_new_totals_dataset:
        scaled_list_int = scale_list_int(new_total, original_list)
        my_dict[new_total] = scaled_list_int
        print("{}: check sum {} : {}".format(new_total, sum(scaled_list_int),scaled_list_int))

    write_to_csv_dict_key_list_values(my_dict, csv_file_name)

    



if __name__ == '__main__':
    main()