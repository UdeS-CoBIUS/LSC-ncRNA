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
    string test_name = "test"; // -tn
    int nb_families = 10; // -nf
    int min_nb_seqs_allowed = 4; // -mins
    int max_nb_seqs_allowed = 1000; // -maxs
    int min_length_motif = 2; // -minl
    int max_length_motif = 10; // -maxl
    int is_delete_subMotifs = 0; // -d (false, true)
    int Beta = 40; // -b
    int Alpha = -1; // -a   // -1, mean don't use the Alpha paramter ==> whatever the variance we accepte it.
    int nbOccrs_allowed = 1; // -g (gamma) , lower bound, 1 is defult

};

void print_args_definition() {
    cout << "-in : <string> a path directory for fasta files" << endl;
    cout << "-tn : <string> test name, give a specifc name to you your exepremnt" << endl;
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
    cout << "test_name: " << arg.test_name << endl;
    cout << "nb_families: " << arg.nb_families << endl;
    cout << "min_nb_seqs_allowed: " << arg.min_nb_seqs_allowed << endl;
    cout << "max_nb_seqs_allowed: " << arg.max_nb_seqs_allowed << endl;
    cout << "min_length_motif: " << arg.min_length_motif << endl;
    cout << "max_length_motif: " << arg.max_length_motif << endl;
    cout << "is_delete_subMotifs: " << arg.is_delete_subMotifs << (is_delete_subMotifs ? " True " : " False ") << endl;
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
        } else if (arg == "-tn") {
            res.test_name = argv[++i];
        } else {
            cout << "Usage: " << argv[0] << " -in <string> -nf <integer> -mins <integer> -maxs <integer> -minl <integer> -maxl <integer> -d <integer> -b <integer> -a <integer> -g <integer>" << endl;
            print_args_definition();
            exit(1);
        }
    }

    return res;
}


int main(int argc, char *argv[]) {

    Args args = get_args(argc, argv);


    // the output matrix in the csv file 

    string output_csv_file = "del_No_nbF_";

    if(args.is_delete_subMotifs) output_csv_file = "del_Yes_nbF_";

    output_csv_file += args.test_name;
    output_csv_file += "_min_" + util::to_string(args.min_length_motif);
    output_csv_file += "_max_" + util::to_string(args.max_length_motif);
    output_csv_file += "_beta_" + util::to_string(args.Beta);
    output_csv_file += "_alpha_" + util::to_string(args.Alpha);
    output_csv_file += "_nbOccrs_" + util::to_string(args.nbOccrs_allowed);


    print_args(args);    

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

    //FastaFilesReader::groupALLFamiliesSeqsInOneFile(dir_name);


    // this construction for no max no min (default min is 1, max is the length of the sequnces)
    SuffixTree_QuadraticTime gst(args.max_length_motif, args.min_length_motif);

    gst.GenerateGeneralizedSuffixTree(list_families_sequences);

    //PrintTree::PrintSuffixTree(gst.getRootSuffixTree());

    CommonMotifs cms(gst, Beta, Alpha, nbOccrs_allowed);


    if (args.is_delete_subMotifs) {
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