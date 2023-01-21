#include <iostream>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <set>
#include <map>
#include <forward_list>
#include "SuffixTree_QuadraticTime.h"
#include "PrintTree.h"
#include "Old_CommonMotifs.h"
#include "FastaFilesReader.h"
#include "CommonMotifs.h"
#include "hasher_Container.h"

#include <unistd.h>
#include <iostream>
using namespace std;

// Struct to store the parsed command-line arguments
struct Args {

    string dir_name;  // -in
    int nb_families = 10; // -nf
    int min_nb_seqs_allowed = 4; // -mins
    int max_nb_seqs_allowed = 1000; // -maxs
    int min_length_motif = 2; // -minl
    int max_length_motif = 10; // -maxl
    int is_delete_subMotifs = 0; // -d (false, true)
    int Beta = 40; // -b
    int Alpha = -1; // -a
    int nbOccrs_allowed = 1; // -g (gamma)

};

void print_args_definition() {
    cout << "-in : <string> a path directory for fasta files" << endl;
    cout << "-nf : <integer> number of families "<< endl;
    cout << "-mins : <integer>, min number of sequences (default 4)"<< endl;
    cout << "-maxs : <integer>, max number of sequences"<< endl;
    cout << "-minl : <integer>, min length of motif" << endl;
    cout << "-maxl : <integer>, max length of motif" << endl;
    cout << "-d : <integer> (0: false, 1 or other: true), is delete sub-motifs" << endl;
    cout << "-b : <integer> beta (between [0 and 100])" << endl;
    cout << "-a : <integer>, alpha  (-1 default no alpha, or: 0 equal number of occurrences, or 1,2,3,.... )" << endl;
    cout << "-g : <integer> ( >=1), gamma, number of occurrences allowed" << endl;
}


void print_args(Args arg) {
    cout << "dir_name: " << arg.dir_name << endl;
    cout << "nb_families: " << arg.nb_families << endl;
    cout << "min_nb_seqs_allowed: " << arg.min_nb_seqs_allowed << endl;
    cout << "max_nb_seqs_allowed: " << arg.max_nb_seqs_allowed << endl;
    cout << "min_length_motif: " << arg.min_length_motif << endl;
    cout << "max_length_motif: " << arg.max_length_motif << endl;
    cout << "is_delete_subMotifs: " << arg.is_delete_subMotifs << endl;
    cout << "Beta: " << arg.Beta << endl;
    cout << "Alpha: " << arg.Alpha << endl;
    cout << "nbOccrs_allowed: " << arg.nbOccrs_allowed << endl;
}

// Function to parse command-line arguments using getopt
Args get_args(int argc, char* argv[]) {

    Args res;
 
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-in") {
            res.dir_name = argv[++i];
        } else if (arg == "-nf") {
            res.nb_families = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-mins") {
            res.min_nb_seqs_allowed = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-maxs") {
            res.max_nb_seqs_allowed = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-minl") {
            res.min_length_motif = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-maxl") {
            res.max_length_motif = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-d") {
            res.is_delete_subMotifs = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-b") {
            res.Beta = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-a") {
            res.Alpha = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-g") {
            res.nbOccrs_allowed = strtol(argv[++i], nullptr, 10);
        } else {
            cout << "Usage: " << argv[0] << " -in <string> -nf <integer> -mins <integer> -maxs <integer> -minl <integer> -maxl <integer> -d <integer> -b <integer> -a <integer> -g <integer>" << endl;
            print_args_definition();
            exit(1);
        }
    }

    return res;
}


int main(int argc, char *argv[]) {

    string dir_name = R"(/data/chei2402/ibra/test_infernal/nbF_all_nbSeqs_min_3/Train)";
    dir_name=argv[1];

    string output_csv_file;
    //string output_csv_file = "nbF_380_nbSeqs_min_2_max_2";
    //string output_csv_file = argv[2];

    uint32_t nb_families = 10;
    uint32_t min_nb_seqs_allowed = 20;
    uint32_t max_nb_seqs_allowed = 40;

    //int nb = 16;
    //for (int nb = 2; nb <= 20; nb+=1)
    //{

    output_csv_file = "del_No_nbF_";

    uint32_t min_length_motif = 1;
    min_length_motif = strtol(argv[2], nullptr, 10);
    uint32_t max_length_motif = 7;
    max_length_motif = strtol(argv[3], nullptr, 10);

    uint32_t Beta = 0;
    Beta = strtol(argv[4], nullptr, 10);

    int Alpha = -1; // -1, mean don't use the Alpha paramter ==> whatever the variance we accepte it.

    unsigned int nbOccrs_allowed = 1; // lower bound allowed, 1 is default
    nbOccrs_allowed = strtol(argv[5], nullptr, 10);

    string families_name = argv[6];

    // string save_csv_mode = argv[7]; // 1, 2, or 3 // according to test settings results. we will use cms.saveMatrixCMS_ToCsv_File_dircrly(output_csv_file);

    uint32_t is_delete_subMotifs = 0;
    //is_delete_subMotifs = strtol(argv[1], nullptr, 10);
    if(is_delete_subMotifs) output_csv_file = "del_Yes_nbF_";

    output_csv_file += families_name;
    output_csv_file += "_min_" + util::to_string(min_length_motif);
    output_csv_file += "_max_" + util::to_string(max_length_motif);
    output_csv_file += "_beta_" + util::to_string(Beta);
    output_csv_file += "_alpha_" + util::to_string(Alpha);
    output_csv_file += "_nbOccrs_" + util::to_string(nbOccrs_allowed);

    //getMainArgv(argc, argv, dir_name, nb_families, min_nb_seqs_allowed, max_nb_seqs_allowed, min_length_motif,max_length_motif, is_delete_subMotifs, Beta);

    cout << "dire name : " << dir_name << endl;
    cout << "nb families : " << nb_families << endl;
    cout << "Min nb seqs : " << min_nb_seqs_allowed << endl;
    cout << "Max nb seqs : " << max_nb_seqs_allowed << endl;
    cout << "Min length motif : " << min_length_motif << endl;
    cout << "Max length motif : " << max_length_motif << endl;
    cout << (is_delete_subMotifs ? " True " : " False ") << is_delete_subMotifs << endl;
    cout << "BETA : " << Beta << endl;
    cout << "Alpha : " << Alpha << endl;
    cout << "nbOccrs_allowed : " << nbOccrs_allowed << endl;

    auto start = chrono::high_resolution_clock::now();

    //FastaFilesReader::getSaveInfosRNAFamiliesCSVFile(dir_name);

    //string str="ACTGactg";
    //vector<string> listStrings={"ttgg","ccttgg","ttggcc"}; // sub motif for all subtree
    //vector<string> listStrings={"ACTGactg"}; // sub motif for all subtree

    //auto list_families_files_names = FastaFilesReader::get_First_N_Files(dir_name,nb_families);
    //auto list_families_sequences = FastaFilesReader::getListFamiliesSequences(list_families_files_names);
    //auto list_families_sequences = FastaFilesReader::getListFamiliesSequences_FirstNFiles(dir_name,nb_families,max_nb_seqs_allowed);
    //auto list_families_sequences = FastaFilesReader::getListFamiliesSequences_FirstNFiles_MinMax(dir_name,nb_families,min_nb_seqs_allowed,max_nb_seqs_allowed);
    auto list_families_sequences = FastaFilesReader::getListFamiliesSequences(dir_name);

    //return 0;
//    int i=1;
//    int total=0;
//    for(const auto &list_family:list_families_sequences)
//    {
//        cout<<"family["<<i++<<"] : "<<list_family.size()<<endl;
//        total +=list_family.size();
//    }
//    cout<<"total : "<<total<<endl;

    //FastaFilesReader::groupALLFamiliesSeqsInOneFile(dir_name);

    //return 0;

    //SuffixTree_QuadraticTime gst;  // this construction for no max no min (default min is 1, max is the length of the sequnces)
    SuffixTree_QuadraticTime gst(max_length_motif, min_length_motif);
    //gst.GenerateSuffixTree(str);
    //gst.GenerateGeneralizedSuffixTree(listStrings);
    gst.GenerateGeneralizedSuffixTree(list_families_sequences);

    //PrintTree::PrintSuffixTree(gst.getRootSuffixTree());

    //uint32_t nb_total_motifs;
    //uint32_t nb_total_motifs_after_deletion_sm;
    //-----------------------  Old_CommonMotifs  --------------------------
    /*
    Old_CommonMotifs cms(gst, Beta);

    auto start = chrono::high_resolution_clock::now();
    cms.generateMatrixCmsSeqs();
    //cms.print_list_cm_umap_seqId_nbOccs();

    nb_total_motifs = cms.list_of_all_commonMotifs.size();

    //cout<<"--------------------"<<endl;
    if(is_delete_subMotifs)
        cms.deletSubMotifsFromTheMatrix();

    nb_total_motifs_after_deletion_sm = cms.list_of_all_commonMotifs.size();

    auto end = chrono::high_resolution_clock::now();
    //cms.print_matrix_commonMotifs_seqsId();

    */

    //-----------------------  New CommonMotifs  --------------------------


    CommonMotifs cms(gst, Beta, Alpha, nbOccrs_allowed);

    ///*
    //auto start = chrono::high_resolution_clock::now();
    //cms.cmsExtractionSelection();
    //cms.print_list_cm_umap_seqId_nbOccs();

    //nb_total_motifs = cms.get_nb_totals_cms();

    //cout<<"--------------------"<<endl;
    //if(is_delete_subMotifs)
    //cms.deleteSubMotifs_at_end_ES();
    //cms.deleteSubMotifs_at_end_ES_usingHash();

    //cms.generateMatrixCmsSeqs();
    //cms.print_matrix_cms_seqsId();

    //auto end = chrono::high_resolution_clock::now();

    //nb_total_motifs_after_deletion_sm = cms.get_nb_totals_cms();

    //*/

    /*

    cms.cmsExtractionSelectionDeletionCSMs();
    cms.generateMatrixCmsSeqs();
    //cms.print_list_cm_umap_seqId_nbOccs();


    //nb_total_motifs_after_deletion_sm = cms.get_nb_totals_cms();
    //cms.print_matrix_cms_seqsId();

     */

    //cms.saveToCsvMatrixCommonMotifsSeqsIds_buffer_reserve(); // 18.71



    if (is_delete_subMotifs) {
        cms.cmsExtractionSelectionDeletionCSMs();
    } else {
        cms.cmsExtractionSelection();
    }

    cms.generateMatrixCmsSeqs();

    auto end_before_saving = chrono::high_resolution_clock::now();

    cms.print_infos();

    // according to test results, we will use saveMatrixCMS_ToCsv_File_dircrly
    /// cms.saveMatrixCMS_ToCsv_File_dircrly(output_csv_file);
    cms.saveMatrixCMS_ToCsv_File(output_csv_file);
    /*
    if(save_csv_mode == "0"){ // default
        cms.saveMatrixCMS_ToCsv_File(output_csv_file);
    }
    else if(save_csv_mode == "0std"){ //
        cms.saveMatrixCMS_ToCsv_File_stdtostring(output_csv_file);
    }
    else if(save_csv_mode == "0stdsize"){ //
        cms.saveMatrixCMS_ToCsv_File_stdtostring_size_estimated(output_csv_file);
    }
    else if(save_csv_mode == "1"){ //
        cms.saveMatrixCMS_ToCsv_File_ofstream_write(output_csv_file);
    }
    else if(save_csv_mode == "2"){ //
        cms.saveMatrixCMS_ToCsv_File_dircrly(output_csv_file);
    }
    else if(save_csv_mode == "3"){ //
        cms.saveMatrixCMS_ToCsv_File_rowByrow(output_csv_file);
    }
    else if(save_csv_mode == "4"){ //
        unsigned int chunk_size = 1000;
        cms.saveMatrixCMS_ToCsv_File_Chunked(output_csv_file, chunk_size);
    }
    */

    auto end = chrono::high_resolution_clock::now();

    //cout<<" All extracted nb_motifs = "<<nb_total_motifs<<endl;
    //cout<<" Remaining nb_motifs After deletion= "<<nb_total_motifs_after_deletion_sm<<endl;

    // Calculating total time taken by the program.

    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end_before_saving - start).count();
    time_taken *= 1e-9;
    cout << "Time taken before save_to_csv is : " << fixed << time_taken << setprecision(9) << " sec" << endl;

    time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Time taken by whole program is : " << fixed << time_taken << setprecision(9) << " sec" << endl;



    // this is special delet a cause de loop
    list_families_sequences.clear();

    //} // end loop for to generate all matrix at the same time.

    return 0;
}