//
// Created by ibra on 12/26/2019.
//

#ifndef SUFFIX_TREE_FASTAFILESREADER_H
#define SUFFIX_TREE_FASTAFILESREADER_H


#include <vector>
#include <string>
//#include <filesystem>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

class FastaFilesReader {

    static const char newline;

public:
    static vector<string> getListSequences(const string &file_name);
    static string getNbMinMaxMeanFamilySequences(const string &file_name);
    static vector<string> getInfosListFamiliesSequences(const string &path_dir);
    static void getSaveInfosRNAFamiliesCSVFile(const string &path_dir);

    static vector<vector<string>> getListFamiliesSequences(const vector<string> &list_families_files_names);

    static vector<vector<string>> getListFamiliesSequences_FirstNFiles(const string &path_dir, uint32_t nb_files, uint32_t nb_seqs_allowed);

    static vector<vector<string>> getListFamiliesSequences(const string &dir_families_files_names);

    static vector<string> getListFiles(const string &path_dir);

    static inline char separator();

    static vector<string> get_Random_N_Files(const string &path_dir, unsigned int nb_files);

    static vector<string> get_First_N_Files(const string &path_dir, unsigned int nb_files);

    static bool isContainLessThanNSeqs(const string &file_name, uint32_t nb_seqs_allowed);

    static bool isNbSeqsBetween_MinMax(const string &file_name, uint32_t min, uint32_t max);

    static vector<vector<string>>
    getListFamiliesSequences_FirstNFiles_MinMax(const string &path_dir, uint32_t nb_files,
                                                uint32_t nb_seqs_min, uint32_t nb_seqs_max);

    static void
    get_Save_N_Random_Family_nbSeqs_in_MinMax(const string &path_dir_in, const string &path_dir_out, uint32_t nb_files,
                                                uint32_t nb_seqs_min, uint32_t nb_seqs_max);

    static void groupALLFamiliesSeqsInOneFile(const string &path_dir);

    static void construct_Train_Test_files(const string &path_dir_in, const string &path_dir_out, uint32_t nb_files,
                                           uint32_t nb_seqs_min, uint32_t nb_seqs_max,
                                           uint32_t percentage_nb_seqs_test);

    static void construct_Train_Test_files(const string &path_dir_in, const string &path_dir_out,
                                           uint32_t nb_seqs_min,
                                           uint32_t percentage_nb_seqs_test);

    static void construct_Train_Test_files(const string &path_dir_in, const string &path_dir_out, uint32_t percentage_nb_seqs_test);

    static void copy_Train_Test_files(const string &path_dir_in, const string &path_dir_out, uint32_t nb_files,
                                           uint32_t nb_seqs_min, uint32_t nb_seqs_max);

    static void write_to_file(const string &buffer, const string &file_name, const string &dir);

    static void creat_dir_c(char *path_dir);

    static string getFileName(string filePath, bool withExtension, char seperator);

    static int getNbSeqs(const string &file_name);

    static void
    get_Families_files(const string &path_dir_in, const string &path_dir_out, uint32_t max_lenght_seq,
                       uint32_t nb_seqs_min,
                       uint32_t nb_seqs_max);

    static bool areAllSeqsLessEqualMax(const string &file_name, uint32_t max);

    static void copyFile(const string &fileNameFrom, const string &fileNameTo);
};


#endif //SUFFIX_Tstatic REE_FASTAFILESREADER_H
