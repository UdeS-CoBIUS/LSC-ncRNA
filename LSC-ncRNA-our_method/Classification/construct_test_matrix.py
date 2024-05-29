import csv
import pandas as pd
import datatable as dt
import ahocorasick
from suffix_trees import STree
import random
import time

import model

import k_mer_Based as kmr

def generate_random_string(str_chars,n):
    return ''.join(random.choices(str_chars, k=n))

def get_matrix_nbOcrrs_listStr_AhoCorasick(list_motifs_in, list_str):
    matrix = []
    list_nboccrs = [0]*len(list_motifs_in)

    trie = ahocorasick.Automaton()
    for idx, m in enumerate(list_motifs_in):
        trie.add_word(m,idx)

    trie.make_automaton()

    for str in list_str:
        for end_idx, key in trie.iter(str):
            list_nboccrs[key] += 1
        matrix.append(list_nboccrs.copy())
        list_nboccrs = [0]*len(list_motifs_in)
        #list_nboccrs = np.zeros(len(list_motifs_in),dtype=int)

    return matrix



def get_list_cm_nbOcrrs_listStr(list_motifs_in, list_str):

    matrix = []
    list_nboccrs = []

    for str in list_str:
        st = STree.STree(str)
        for cm in list_motifs_in:
            list_nboccrs.append(len(st.find_all(cm)))
        matrix.append(list_nboccrs.copy())
        list_nboccrs.clear()

    #print(" Matrix[5] : ------------------------------")
    #print(matrix[5])
    return matrix


def get_list_cm_nbOcrrs(list_motifs_in, str):

    list_cm_nboccrs = []

    st = STree.STree(str)

    for cm in list_motifs_in:
        list_cm_nboccrs.append(len(st.find_all(cm)))

    return list_cm_nboccrs

def test_SuffixTree_vs_AhoCorsick():
    str_RNA = "ACUG"
    size_str = 250
    nb_str = 1000
    list_strs = [None] * nb_str
    for i in range(nb_str):
        list_strs[i] = generate_random_string(str_RNA, size_str)

    print(list_strs[1])

    print(" size of string = ", size_str)
    print("nb strings = ", len(list_strs))

    list_motifs = kmr.generate_kMer_MinMax(max=8)

    print(" nb motifs = ", len(list_motifs))

    print(" wait for AhoCorasick....")
    start_time = time.time()
    matrix_AC = get_matrix_nbOcrrs_listStr_AhoCorasick(list_motifs, list_strs)
    end_time = time.time()

    print(" time for AC = ", end_time - start_time)

    print(" wait for Suffixe Tree....")
    start_time = time.time()
    matric_ST = get_list_cm_nbOcrrs_listStr(list_motifs, list_strs)
    end_time = time.time()

    print(" time for ST = ", end_time - start_time)



    print(" is tow matriw same : ", matric_ST == matrix_AC)


def check_csvMatrix_genrateMatrix_Train(csv_file_path, dir_in_train_files, file_ext):
    print(" in check_csvMatrix_genrateMatrix_Train ...")

    start_time = time.time()
    data_arn = pd.read_csv(csv_file_path, engine='c') # On importe la dataset
    #names_ = [str]*80822
    #dtype_ = [int]*80822
    #data_arn = pd.read_csv(csv_file_path , engine='c', names=str,dtype=int) # On importe la dataset
    #data_arn = pd.read_csv(csv_file_path , engine='c', names=names_,dtype=dtype_) # On importe la dataset
    end_time = time.time()
    print(" time Read CSV File = ",end_time-start_time, " s")
    print("pd.read_csv ...")
    data_arn = data_arn.rename(columns={"familyId": "Classe"})
    data_arn = data_arn.drop(['seqIdInFam', 'index'], axis=1)
    data_arn.drop(data_arn.filter(regex="Unname"), axis=1, inplace=True)
    data_arn = data_arn.astype(int) # On change le type to int
    print(" rename, drop, filter, change type...")

    data_classe = data_arn["Classe"]
    data_motif = data_arn.drop(['Classe'], axis=1)
    list_motifs = list(data_motif) # save list of use motifs in classification

    print(" len( Data_motifs ) = ",len(data_motif))
    list_all_seqs = model.Model.get_all_seqs(dir_in_train_files, file_ext)
    list_all_classes = model.Model.get_all_classes(dir_in_train_files, file_ext)

    print("list all seqs len: ", len(list_all_seqs))
    print(" is len all seqs = len all calsse : ", len(list_all_seqs) == len(list_all_classes))

    generated_matrix = get_matrix_nbOcrrs_listStr_AhoCorasick(list_motifs, list_all_seqs)

    dftest = pd.DataFrame(generated_matrix, columns=list_motifs)

    print(" is matrix CSV == Generated from Seq", data_motif == dftest)

    isSame = True
    for row in dftest:
        if row not in data_motif:
            isSame = False
            break
    print(" isSame = ", isSame)


def load_csv_manualy(path_csv_file):
    #csv.field_size_limit(2_147_483_647)
    csv.field_size_limit(2147483647)
    nb_row = 0
    matrix = []
    with open(path_csv_file, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for idx_r,row in enumerate(spamreader):
            nb_row += 1

            if idx_r == 0:
                matrix.append(row)
            else:
                vect = [int] * len(row)
                for idx,i in enumerate(row):
                    vect[idx] = int(i)
                matrix.append(vect)
    print(" nb row = ", nb_row)
    return matrix


def test_load_csv():
    #test_csv_matrix = R"C:\Users\ibra\OneDrive - USherbrooke\Project\MotifsExtractionSelection\cmake-build-debug\test_csv_file_matrix_50_Copy.csv"
    test_csv_matrix = R"C:\Users\ibra\Desktop\Infernal\test_200F.csv"

    start_time = time.time()
    #matrix_p=load_csv_manualy(test_csv_matrix)
    #matrix = pd.DataFrame(matrix_p)
    #matrix = pd.read_csv(test_csv_matrix)
    #matrix = pd.read_csv(test_csv_matrix,engine='c')
    #matrix = pd.read_csv(test_csv_matrix,engine='c',sep=',')
    #matrix = pd.read_csv(test_csv_matrix,engine='c',sep=',', low_memory=False)

    #matrix = pd.read_csv(test_csv_matrix,engine='c',sep=',',dtype=int)
    #TextFileReaderObject = pd.read_csv(test_csv_matrix,engine='c',sep=',',dtype=int,chunksize=100_000)
    #matrix = pd.concat(chunk for chunk in TextFileReaderObject)

    #dask_dataframe = dask.dataframe.read_csv(test_csv_matrix)
    #matrix = dask_dataframe.compute()

    data = dt.fread(test_csv_matrix)
    #matrix = data.to_pandas().astype('int32')
    matrix = data.to_pandas()

    #print(matrix.iloc[16:26,-10:])
    end_time = time.time()
    print(" time = ", end_time - start_time)
    #print(" nb row = ", len(matrix))
    #print(" nb row = ", len(data))
    print(data.head(5))
    print(matrix.head())
