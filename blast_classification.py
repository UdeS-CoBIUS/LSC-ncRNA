import math
import sys
import time
import subprocess
from Bio.Blast import NCBIXML


def create_blast_db(in_files_dir, name_db):
    # 1) concatenate all file in one file
    fasta_file_out = "fasta_files_" + name_db
    cat_cmd = "cat " + "\'" + in_files_dir + "\'" + "* > " + fasta_file_out
    print(cat_cmd)
    subprocess.call(cat_cmd, shell=True)

    # 2) create the database
    db_cmd = "makeblastdb  -in " + fasta_file_out + " -dbtype nucl -title " + name_db + " -out " + name_db
    print(db_cmd)
    subprocess.call(db_cmd, shell=True)


def search_blast_db(test_in_file_dir, name_db):
    # 1) concatenate all file in one file
    query_input_file = "query_" + name_db
    cat_cmd = "cat " + "\'" + test_in_file_dir + "\'" + "* > " + query_input_file
    # print(cat_cmd)
    subprocess.call(cat_cmd, shell=True)

    # using -outfmt 5 : pour get xml result file.
    # -word_size int_value : the min size of the word search by blast.

    # A) by default
    # search_db_cmd = "blastn -db " + name_db + " -query " + query_input_file + " -out res_" + query_input_file + ".xml -outfmt 5"

    # B) change word_size
    search_db_cmd = "blastn -db " + name_db + " -query " + query_input_file + " -out res_wz7_" + query_input_file + ".xml -outfmt 5 -word_size 7"
    # print(search_db_cmd)

    subprocess.call(search_db_cmd, shell=True)


def get_train_path(tm):
    if tm == "rej_all":
        return R"/data/chei2402/ibra/test_infernal/deep_ncrna_datasets/original_family_rejection/Train/"

    if tm == "rej_rd":
        return R"/data/chei2402/ibra/test_infernal/deep_ncrna_datasets/original_family_rejection_random/Train/"

    if tm == "rej_sh":
        return R"/data/chei2402/ibra/test_infernal/deep_ncrna_datasets/original_family_rejection_shuffle/Train/"

    return R"/data/chei2402/ibra/test_infernal/deep_ncrna_datasets/original_family_rejection_dinucliutide/Train/"


def get_test_path(tm):
    if tm == "rej_all":
        return R"/data/chei2402/ibra/test_infernal/deep_ncrna_datasets/original_family_rejection/Test/"

    if tm == "rej_rd":
        return R"/data/chei2402/ibra/test_infernal/deep_ncrna_datasets/original_family_rejection_random/Test/"

    if tm == "rej_sh":
        return R"/data/chei2402/ibra/test_infernal/deep_ncrna_datasets/original_family_rejection_shuffle/Test/"

    return R"/data/chei2402/ibra/test_infernal/deep_ncrna_datasets/original_family_rejection_dinucliutide/Test/"


def create_blast_db_for_multiple_benchmark(list_size_families):
    # 1) creating blast local data base
    for i in list_size_families:
        train_in_file_dir = get_train_path(i)
        name_db = "db_deep_ncrna_bn_" + str(i) + "_train"

        start = time.time()
        create_blast_db(train_in_file_dir, name_db)
        end = time.time()
        print(" dir {} Creating in : {} sec ..................".format(i, end - start))


def search_blast_db_for_multiple_benchmark(list_size_families):
    # 2) searching the blast local data base
    for i in list_size_families:
        test_in_file_dir = get_test_path(i)
        name_db = "db_deep_ncrna_bn_" + str(i) + "_train"

        start = time.time()
        search_blast_db(test_in_file_dir, name_db)
        end = time.time()
        print(" dir {} Searching in : {} sec ..................".format(i, end - start))


def remove_extension(in_str):
    pos = in_str.rfind(".")  # get pos of the last "."
    if pos != -1:
        return in_str[:pos]
    return in_str


# 1) NNCM : Naif Naif Calss Method =>
# seq query goes with family which is nb of seqs majoritaire , donc > total nb seqs /2
# dans notre cas, on test seulement cette condition pour notre query , si elle est bien classer ou pas
# nb_res_same_family : calculer le nombre des seqs res de la meme family que seq query
# si (nb_res_same_family) est majoritaire, donc  il est > (total number /2), on dit que la prediction est correct
def parse_result_classification_naif_naif(blast_xml_result):
    result_handle = open(blast_xml_result)
    blast_records = NCBIXML.parse(result_handle)
    # blast_records = list(NCBIXML.parse(result_handle))

    nb_true_classified = 0
    total_nb_query = 0

    for rec in blast_records:
        id_query = remove_extension(rec.query)
        same = 0

        for ali in rec.alignments:
            if id_query == remove_extension(ali.hit_def):
                same = same + 1

        if len(rec.alignments) / 2 < same:
            nb_true_classified = nb_true_classified + 1

        total_nb_query = total_nb_query + 1

    print("total_nb_query = ", total_nb_query)
    print("nb_true_classified = ", nb_true_classified)
    print(" accuracy = ", nb_true_classified / total_nb_query)


# 2) NCM : Naif Naif Calss Method =>
# Cette method est calairement pas bonne, elle marche seuelet si on a deux type de seqs dans deux famillies, donc avoir plus que la moitie est logique, mais si on plusieurs? Sa ne marche pas, on doit trouver quele famillie est majoritaire
# s'est simple, compter le nb de seqs de chaque famille, et prend le plus grand.
def parse_result_classification_naif(blast_xml_result):
    result_handle = open(blast_xml_result)
    blast_records = NCBIXML.parse(result_handle)
    # blast_records = list(NCBIXML.parse(result_handle))

    list_query = []

    nb_true_classified = 0
    total_nb_query = 0
    total_nb_families_for_all_queries = 0

    for rec in blast_records:

        id_query = remove_extension(rec.query)
        dict_query_found_seq = {}
        # count nb seqs by family
        for ali in rec.alignments:
            res_seq_id = remove_extension(ali.hit_def)
            if res_seq_id not in dict_query_found_seq:
                dict_query_found_seq[res_seq_id] = 1
            else:
                dict_query_found_seq[res_seq_id] += 1

        # debug:
        total_nb_families_for_all_queries += len(dict_query_found_seq)
        # get family with highest nb of seqs
        max = 0
        best_fam = None
        for family, nb_seq in dict_query_found_seq.items():
            if nb_seq > max:
                max = nb_seq
                best_fam = family
        # here , best_fam is the classification for id_query. because we test, if we get a good classification, we test their id after

        # nb true classification
        if id_query == best_fam:
            nb_true_classified = nb_true_classified + 1

        # count total query
        total_nb_query = total_nb_query + 1

    print("total_nb_query = ", total_nb_query)
    print("nb_true_classified = ", nb_true_classified)
    print(" accuracy = ", nb_true_classified / total_nb_query)
    print(" Average nb families by query : ", total_nb_families_for_all_queries / total_nb_query)


# 3) NCM_Score : Naif Naif Calss Method + score => sum all socres for each family, the family with max score win.
# 3) NCM + Score
# Using method NCM, if f1 have 9 seq , and f2 have 10 seqs, in this case the calssification goes for f2
# but,  if seq of f1 have better score than of seqs of f2
# So, we have to include score + nb seqs ,  we have to come up with an eqution for that
# naif score: we can just relay on average score to classify. Exemple f1 with 9 seqs, avrage score is best than f2 with 20 seqs. So classification goes for f1.  but f2 have 20 seqs, and this have a sinification
# nb seqs + score : we can sum all the scores from each seqs,  and the family that have max score win.

def parse_result_classification_NCM_Score(blast_xml_result):
    result_handle = open(blast_xml_result)
    blast_records = NCBIXML.parse(result_handle)
    # blast_records = list(NCBIXML.parse(result_handle))

    list_query = []

    nb_true_classified = 0
    total_nb_query = 0
    total_nb_families_for_all_queries = 0

    for rec in blast_records:

        id_query = remove_extension(rec.query)
        dict_query_found_seq = {}
        # count nb seqs by family
        for ali in rec.alignments:
            res_seq_id = remove_extension(ali.hit_def)
            score = ali.hsps[0].score
            if res_seq_id not in dict_query_found_seq:
                # dict_query_found_seq[res_seq_id] = 1
                dict_query_found_seq[res_seq_id] = score
            else:
                # dict_query_found_seq[res_seq_id] += 1
                dict_query_found_seq[res_seq_id] += score

        # debug:
        total_nb_families_for_all_queries += len(dict_query_found_seq)
        # get family with highest nb of seqs
        max = 0
        best_fam = None
        # for family, nb_seq in dict_query_found_seq.items():
        for family, family_score in dict_query_found_seq.items():
            if family_score > max:
                max = family_score
                best_fam = family
        # here , best_fam is the classification for id_query. because we test, if we get a good classification, we test their id after

        # nb true classification
        if id_query == best_fam:
            nb_true_classified = nb_true_classified + 1

        # count total query
        total_nb_query = total_nb_query + 1

    print("total_nb_query = ", total_nb_query)
    print("nb_true_classified = ", nb_true_classified)
    print(" accuracy = ", nb_true_classified / total_nb_query)
    print(" Average nb families by query : ", total_nb_families_for_all_queries / total_nb_query)


# this is the same of : parse_result_classification_NCM_Score
# but we compute the average
# I espect that work not god as the only sum of score,
# because, it doesn't include the information of nb seqs by family
# example :
# f1 with 10 seqs and average score of 10 (sum = 100)
# f2 with 20 seqs and average score of 9 (sum = 180)
# in previuos :parse_result_classification_NCM_Score, f2 win
# in this, f1 win

def parse_result_classification_NCM_Score_average(blast_xml_result):
    start_time = time.time()
    result_handle = open(blast_xml_result)
    blast_records = NCBIXML.parse(result_handle)
    # blast_records = list(NCBIXML.parse(result_handle))

    list_query = []

    nb_true_classified = 0
    total_nb_query = 0
    total_nb_families_for_all_queries = 0

    for rec in blast_records:

        id_query = remove_extension(rec.query)
        dict_query_found_seq = {}
        # count nb seqs by family
        for ali in rec.alignments:
            res_seq_id = remove_extension(ali.hit_def)
            score = ali.hsps[0].score
            if res_seq_id not in dict_query_found_seq:
                # dict_query_found_seq[res_seq_id] = 1
                dict_query_found_seq[res_seq_id] = [score]
            else:
                # dict_query_found_seq[res_seq_id] += 1
                dict_query_found_seq[res_seq_id].append(score)

        # debug:
        total_nb_families_for_all_queries += len(dict_query_found_seq)
        # get family with highest nb of seqs
        max = 0
        best_fam = None
        for family, list_family_score in dict_query_found_seq.items():
            average_score = sum(list_family_score) / len(list_family_score)
            if average_score > max:
                max = average_score
                best_fam = family
        # here , best_fam is the classification for id_query. because we test, if we get a good classification, we test their id after

        # nb true classification
        if id_query == best_fam:
            nb_true_classified = nb_true_classified + 1

        # count total query
        total_nb_query = total_nb_query + 1

    end_time = time.time()
    print("-------- Result parsing : ", blast_xml_result)
    print("total_nb_query = ", total_nb_query)
    print("nb_true_classified = ", nb_true_classified)
    print(" accuracy = ", nb_true_classified / total_nb_query)
    print(" Average nb families by query : ", total_nb_families_for_all_queries / total_nb_query)
    print(" time : ", (end_time - start_time), " sec")


def parse_result_classification_NCM_Evalue(blast_xml_result):
    start_time = time.time()

    result_handle = open(blast_xml_result)
    blast_records = NCBIXML.parse(result_handle)
    # blast_records = list(NCBIXML.parse(result_handle))

    list_query = []

    nb_true_classified = 0
    total_nb_query = 0
    total_nb_families_for_all_queries = 0

    for rec in blast_records:

        id_query = remove_extension(rec.query)
        dict_query_found_seq = {}
        # count nb seqs by family
        for ali in rec.alignments:
            res_seq_id = remove_extension(ali.hit_def)
            expect = ali.hsps[0].expect
            if res_seq_id not in dict_query_found_seq:
                # dict_query_found_seq[res_seq_id] = 1
                dict_query_found_seq[res_seq_id] = expect
            else:
                # dict_query_found_seq[res_seq_id] += 1
                dict_query_found_seq[res_seq_id] += expect

        # debug:
        total_nb_families_for_all_queries += len(dict_query_found_seq)
        # get family with highest nb of seqs
        max_val = sys.maxsize
        best_fam = None
        # for family, nb_seq in dict_query_found_seq.items():
        for family, family_expect in dict_query_found_seq.items():
            if family_expect < max_val:
                max_val = family_expect
                best_fam = family
        # here , best_fam is the classification for id_query. because we test, if we get a good classification, we test their id after

        # nb true classification
        if id_query == best_fam:
            nb_true_classified = nb_true_classified + 1

        # count total query
        total_nb_query = total_nb_query + 1

    end_time = time.time()
    print("total_nb_query = ", total_nb_query)
    print("nb_true_classified = ", nb_true_classified)
    print(" accuracy = ", nb_true_classified / total_nb_query)
    print(" Average nb families by query : ", total_nb_families_for_all_queries / total_nb_query)
    print(" time : {} sec".format(end_time - start_time))


def parse_result_classification_NCM_Evalue_average(blast_xml_result):
    start_time = time.time()
    result_handle = open(blast_xml_result)
    blast_records = NCBIXML.parse(result_handle)
    # blast_records = list(NCBIXML.parse(result_handle))

    list_query = []

    nb_true_classified = 0
    total_nb_query = 0
    total_nb_families_for_all_queries = 0

    for rec in blast_records:

        id_query = remove_extension(rec.query)
        dict_query_found_seq = {}
        # count nb seqs by family
        for ali in rec.alignments:
            res_seq_id = remove_extension(ali.hit_def)
            expect = ali.hsps[0].expect
            if res_seq_id not in dict_query_found_seq:
                # dict_query_found_seq[res_seq_id] = 1
                dict_query_found_seq[res_seq_id] = [expect]
            else:
                # dict_query_found_seq[res_seq_id] += 1
                dict_query_found_seq[res_seq_id].append(expect)

        # debug:
        total_nb_families_for_all_queries += len(dict_query_found_seq)
        # get family with highest nb of seqs
        max_val = sys.maxsize
        best_fam = None
        for family, list_family_expect in dict_query_found_seq.items():
            average_expect = sum(list_family_expect) / len(list_family_expect)
            if average_expect < max_val:
                max_val = average_expect
                best_fam = family
        # here , best_fam is the classification for id_query. because we test, if we get a good classification, we test their id after

        # nb true classification
        if id_query == best_fam:
            nb_true_classified = nb_true_classified + 1

        # count total query
        total_nb_query = total_nb_query + 1

    end_time = time.time()
    print("total_nb_query = ", total_nb_query)
    print("nb_true_classified = ", nb_true_classified)
    print(" accuracy = ", nb_true_classified / total_nb_query)
    print(" Average nb families by query : ", total_nb_families_for_all_queries / total_nb_query)
    print(" time : ", (end_time - start_time), " sec")


# 4) k nearest neighbors:
# based classification on k best elements in families.
# if tow families, have the same k elements, we go to K+1, on so on.  If they have the same in all sequences, so there are in the same time right classification with same score.
# How to do that:
# 1) get all list of families, with seqs and scores, as we do already
# 2) sort them
# 3) sum the first k element
# 4) take the best.
def parse_result_classification_knn(blast_xml_result, k):
    result_handle = open(blast_xml_result)
    blast_records = NCBIXML.parse(result_handle)
    # blast_records = list(NCBIXML.parse(result_handle))

    list_query = []

    nb_true_classified = 0
    total_nb_query = 0
    total_nb_families_for_all_queries = 0

    for rec in blast_records:

        id_query = remove_extension(rec.query)
        dict_query_found_seq = {}
        # 1) get all list of families, with seqs and scores
        for ali in rec.alignments:
            res_seq_id = remove_extension(ali.hit_def)
            score = ali.hsps[0].score
            if res_seq_id not in dict_query_found_seq:
                dict_query_found_seq[res_seq_id] = [score]
            else:
                dict_query_found_seq[res_seq_id].append(score)

        # debug:
        total_nb_families_for_all_queries += len(dict_query_found_seq)
        # get family with highest nb of seqs
        max_val = 0
        best_fam = None
        for family, family_list_score in dict_query_found_seq.items():

            if len(
                    family_list_score) < k:  # this to avoid bug if len is small than k, also we need exactly k or more elements
                continue

            family_list_score.sort()

            sum = 0
            for i in range(0, k):
                sum += family_list_score[i]

            if sum > max_val:
                max_val = sum
                best_fam = family
        # here , best_fam is the classification for id_query. because we test, if we get a good classification, we test their id after

        # nb true classification
        if id_query == best_fam:
            nb_true_classified = nb_true_classified + 1

        # count total query
        total_nb_query = total_nb_query + 1

    print("total_nb_query = ", total_nb_query)
    print("nb_true_classified = ", nb_true_classified)
    print(" accuracy = ", nb_true_classified / total_nb_query)
    print(" Average nb families by query : ", total_nb_families_for_all_queries / total_nb_query)


# list_score_left and list_score_right : are sorted
# k elements are the already the same, from 0 to k-1
def continue_check_for_best_family(fam_left, list_score_left, fam_right, list_score_right, k):
    min_len = min(len(list_score_left), len(list_score_right))

    for i in range(k, min_len):
        if list_score_left[i] > list_score_right[i]:
            return fam_left
        elif list_score_left[i] < list_score_right[i]:
            return fam_right
        else:
            continue

    # if we arrived to end of one list, without finding a winner, this that the list that have remaining element win
    if len(list_score_left) > min_len:
        return fam_left

    if len(list_score_right) > min_len:
        return fam_right

    # in the case where both len of list == k, we just choose one as winner, the first for example.
    # print(fam_left, len(list_score_left), fam_right, len(list_score_right), " k = ", k)
    return fam_left


# the first take the first best (if there are many that have the same k element best)
# this method, if two have the same k best elements, we go to check the k+1
# ccb : continue check for the best
def parse_result_classification_knn_ccb(blast_xml_result, k):
    result_handle = open(blast_xml_result)
    blast_records = NCBIXML.parse(result_handle)
    # blast_records = list(NCBIXML.parse(result_handle))

    list_query = []

    nb_true_classified = 0
    total_nb_query = 0
    total_nb_families_for_all_queries = 0

    for rec in blast_records:

        id_query = remove_extension(rec.query)
        dict_query_found_seq = {}
        # 1) get all list of families, with seqs and scores
        for ali in rec.alignments:
            res_seq_id = remove_extension(ali.hit_def)
            score = ali.hsps[0].score
            if res_seq_id not in dict_query_found_seq:
                dict_query_found_seq[res_seq_id] = [score]
            else:
                dict_query_found_seq[res_seq_id].append(score)

        # debug:
        total_nb_families_for_all_queries += len(dict_query_found_seq)
        # get family with highest nb of seqs
        max_val = 0
        best_fam = None
        same_best_fam = None
        same_max_val = 0
        for family, family_list_score in dict_query_found_seq.items():

            if len(
                    family_list_score) < k:  # this to avoid bug if len is small than k, also we need exactly k or more elements
                continue

            family_list_score.sort()

            sum_score = 0
            for i in range(0, k):
                sum_score += family_list_score[i]

            if sum_score > max_val:  # get the first best value
                max_val = sum_score
                best_fam = family
            elif sum_score == max_val:  # get the second best value
                same_max_val = sum_score
                same_best_fam = family
            else:
                continue

        # if the first max_val and the second same_max_val, still the same (because it may the first value change)
        # if max_val == same_max_val and max_val is not 0:  # max_val is not 0 : in the case where there is no result, check the best between them ,
        if max_val == same_max_val and max_val != 0:  # max_val is not 0 : in the case where there is no result, check the best between them ,

            # debug:
            # print("best_fam = {},  score = {}".format(best_fam, max_val))
            # print("same_bes = {},  score = {}".format(same_best_fam, same_max_val))

            best_fam = continue_check_for_best_family(best_fam, dict_query_found_seq[best_fam],
                                                      same_best_fam, dict_query_found_seq[same_best_fam], k)

            # print("best_fam = {}".format(best_fam))
            # print("---------------------------------------------")
        # here , best_fam is the classification for id_query. because we test, if we get a good classification, we test their id after

        # nb true classification
        if id_query == best_fam:
            nb_true_classified = nb_true_classified + 1

        # count total query
        total_nb_query = total_nb_query + 1

    print("total_nb_query = ", total_nb_query)
    print("nb_true_classified = ", nb_true_classified)
    print(" accuracy = ", nb_true_classified / total_nb_query)
    print(" Average nb families by query : ", total_nb_families_for_all_queries / total_nb_query)


def parse_classification_multiple(list_size_families):
    for i in list_size_families:
        # name_test_parse = "res_query_db_nbF"+str(i)+"_train.xml"
        name_test_parse = "res_wz7_query_db_deep_ncrna_bn_" + str(i) + "_train.xml"

        # start = time.time()
        # parse_result_classification_naif_naif(name_test_parse)
        # parse_result_classification_naif(name_test_parse)
        # parse_result_classification_NCM_Score(name_test_parse)
        parse_result_classification_NCM_Score_average(name_test_parse)
        # parse_result_classification_NCM_Evalue(name_test_parse)
        # parse_result_classification_NCM_Evalue_average(name_test_parse)
        # for j in [2, 4, 6, 8, 10]:
        #     print(" ----------------- iter : ", j)
        #     start = time.time()
        #     # parse_result_classification_knn(name_test_parse, j)
        #     parse_result_classification_knn_ccb(name_test_parse, j)
        #     end = time.time()
        #     print(" dir {} Searching in : {} sec ..................".format(i, end - start))


def exprement_deep_ncrna_bn():
    mode = sys.argv[1]

    list_size_families = ["rej_all", "rej_rd", "rej_sh", "rej_din"]
    if mode == "-c":
        create_blast_db_for_multiple_benchmark(list_size_families)
    elif mode == "-s":
        search_blast_db_for_multiple_benchmark(list_size_families)
    elif mode == "-p":
        parse_classification_multiple(list_size_families)


def experiment_secondary_and_no_secondary_three_families():
    list_paths_families = [R"/data/chei2402/ibra/test_infernal/small_Lnc_RNAnoSecondary_train_test",
                           R"/data/chei2402/ibra/test_infernal/SecondaryBased_lessSimilarSeqs_train_test",
                           R"/data/chei2402/ibra/test_infernal/Secondary_and_noSecondary"]

    for experiment_path in list_paths_families:
        experiment_secondary_and_no_secondary_one_family(experiment_path)


def experiment_secondary_and_no_secondary_one_family(experiment_path):
    data_set_name = experiment_path.split("/")[-1]  # get the last name of the path

    blast_search_result_name = "res_wz7_query_db_" + data_set_name + "_train.xml"

    train_in_file_dir = experiment_path + "/Train/"
    test_in_file_dir = experiment_path + "/Test/"
    name_db = "db_" + data_set_name + "_train"

    # 1) create blast db
    start = time.time()
    create_blast_db(train_in_file_dir, name_db)
    end = time.time()
    print(" dir {} Creating in : {} sec ..................".format(data_set_name, end - start))

    # 2) search blast db
    start = time.time()
    search_blast_db(test_in_file_dir, name_db)
    end = time.time()
    print(" dir {} Searching in : {} sec ..................".format(data_set_name, end - start))

    # 3) parse classification
    start = time.time()
    parse_result_classification_NCM_Score_average(blast_search_result_name)
    end = time.time()
    print(" dir {} Parsing for classification in : {} sec ..................".format(data_set_name, end - start))


def change_seq_ids():
    in_filepath = R"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_family_train_test_name\Test\CL00001_test.fasta"
    ou_filepath = R"C:\Users\ibra\Desktop\Infernal\Clans ncRNA\Clans_family_train_test_name\Test\CL00001_test_2.fasta"
    new_base_name = "CL00001"
    seq_id = 0
    with open(in_filepath, "r") as file_input:
        with open(ou_filepath, "w") as output:
            for line in file_input:
                if line.startswith('>'):
                    line = ">{}.seq_id_{}\n".format(new_base_name, seq_id)
                    seq_id += 1
                output.write(line)


def experiments(path_train_folder, path_test_folder, data_set_name):
    """
    :param path_train_folder: the path to the directory that contains Train families, each family in separate file.
    :param path_test_folder:  the path to the directory that contains Test  families, each family in separate file.
    :param data_set_name:     the data ste name, or the experiment name
    :return: None
    """

    blast_search_result_name = "res_wz7_query_db_" + data_set_name + "_train.xml"

    name_db = "db_" + data_set_name + "_train"

    # 1) create blast db
    start = time.time()
    create_blast_db(path_train_folder, name_db)
    end = time.time()
    print(" dir {} Creating in : {} sec ..................".format(data_set_name, end - start))

    # 2) search blast db
    start = time.time()
    search_blast_db(path_test_folder, name_db)
    end = time.time()
    print(" dir {} Searching in : {} sec ..................".format(data_set_name, end - start))

    # 3) parse classification
    start = time.time()
    parse_result_classification_NCM_Score_average(blast_search_result_name)
    end = time.time()
    print(" dir {} Parsing for classification in : {} sec ..................".format(data_set_name, end - start))


def main():
    # exprement_deep_ncrna_bn()
    # experiment_secondary_and_no_secondary_three_families()

    path_train_folder = sys.argv[1] # path train dolder: .../Train
    path_test_folder = sys.argv[2] # test folder as : ..../Test
    data_set_name = sys.argv[3] # for example: Clan, secondary, allRfam ... etc

    experiments(path_train_folder, path_test_folder, data_set_name)

    # just for debug
    # blast_xml_result = "res_5.xml"
    # parse_result_classification_NCM_Evalue(blast_xml_result)


def test_code():
    my_list = [1, 2, 3, 4, 8, 5, 6, 7, 8, 0, 9]
    first_max = 0
    idx_1 = None
    same_scd_max = 0
    idx_s = None

    for idx, val in enumerate(my_list):

        if val > first_max:
            first_max = val
            idx_1 = idx
        elif val == first_max:
            same_scd_max = val
            idx_s = idx
        else:
            continue

    print(" first max = ", first_max, " idx = ", idx_1)
    print(" same  max =", same_scd_max, " idx = ", idx_s)

    if first_max == same_scd_max:
        print(" they are the same ")
    else:
        print(" No, ther are not the same ")


if __name__ == '__main__':
    main()
    # test_code()
    #change_seq_ids()
