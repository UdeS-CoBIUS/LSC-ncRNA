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
#include <random> // added for debug
#include <fstream> // added for debug

using namespace std;

// Struct to store the parsed command-line arguments
struct Args {

    string dir_name;  // -in
    string test_name = "test"; // -tn
    int nb_families = -1; // -nf
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
    cout << "is_delete_subMotifs: " << (arg.is_delete_subMotifs ? " True " : " False ") << endl;
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

// generate the name of the output csv file:
string generate_output_csv_file_name(Args args) {
    stringstream ss;
    ss << args.test_name
       << "_nbF_" << args.nb_families
       << "_is_del_" << (args.is_delete_subMotifs ? "yes" : "no")
       << "_min_" << args.min_length_motif
       << "_max_" << args.max_length_motif
       << "_beta_" << args.Beta
       << "_alpha_" << args.Alpha
       << "_nbOccrs_" << args.nbOccrs_allowed
       << ".csv";
    return ss.str();
}

void generate_save_random_matrix(string output_csv_file, int nb_families) {
    // Generate random matrix of 100x100
    const int rows = 10;
    const int cols = nb_families;
    vector<vector<int>> matrix(rows, vector<int>(cols));
    vector<string> items(cols);

    // Random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> len_dist(3, 10);
    uniform_int_distribution<> char_dist(0, 25);
    uniform_int_distribution<> occr_dist(0, 10);

    // Generate random items (substrings)
    for (int j = 0; j < cols; ++j) {
        int len = len_dist(gen);
        string item;
        for (int k = 0; k < len; ++k) {
            item += static_cast<char>('A' + char_dist(gen));
        }
        items[j] = item;
    }

    // Generate random occurrences
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = occr_dist(gen);
        }
    }

    // Save matrix to CSV file
    ofstream outFile(output_csv_file);
    if (!outFile.is_open()) {
        cerr << "Error: Unable to open file " << output_csv_file << endl;
        return;
    }

    // Write header
    outFile << "Index";
    for (const auto& item : items) {
        outFile << "," << item;
    }
    outFile << endl;

    // Write data
    for (int i = 0; i < rows; ++i) {
        outFile << i;
        for (int j = 0; j < cols; ++j) {
            outFile << "," << matrix[i][j];
        }
        outFile << endl;
    }

    outFile.close();
    cout << "Random matrix saved to " << output_csv_file << endl;

    cout<<" Total nb_motifs = "<<cols<<endl;
}

CommonMotifs generate_common_motifs_matrix(Args args){

    std::cout << "we are in generate_common_motifs_matrix" << std::endl;

    //auto list_families_files_names = FastaFilesReader::get_First_N_Files(dir_name,nb_families);
    //auto list_families_sequences = FastaFilesReader::getListFamiliesSequences(list_families_files_names);
    //auto list_families_sequences = FastaFilesReader::getListFamiliesSequences_FirstNFiles(dir_name,nb_families,max_nb_seqs_allowed);
    //auto list_families_sequences = FastaFilesReader::getListFamiliesSequences_FirstNFiles_MinMax(dir_name,nb_families,min_nb_seqs_allowed,max_nb_seqs_allowed);
    auto list_families_sequences = FastaFilesReader::getListFamiliesSequences(args.dir_name);

    for (const auto& list_seq : list_families_sequences) {
        std::cout << "Fam list seq size = " << list_seq.size() << std::endl;
    }

    // this construction for no max no min (default min is 1, max is the length of the sequnces)
    SuffixTree_QuadraticTime gst(args.max_length_motif, args.min_length_motif);

    
    gst.GenerateGeneralizedSuffixTree(list_families_sequences);
    
    // print nb_all_sequences
    std::cout << "nb_all_sequences after gst.GenerateGeneralizedSuffixTree: " << gst.getNbAllSequences() << std::endl;

    // Print the suffix tree
    ///PrintTree::PrintSuffixTree(gst.getRootSuffixTree());

    CommonMotifs cms(gst, args.Beta, args.Alpha, args.nbOccrs_allowed);


    if (args.is_delete_subMotifs) {
        cms.cmsExtractionSelectionDeletionCSMs();
    } else {
        cms.cmsExtractionSelection();
    }

    cms.generateMatrixCmsSeqs();

    // this is special delete, used when we have a loop in the suffix tree
    /// list_families_sequences.clear();
    //} // end loop for to generate all matrix at the same time.

    // debug print nb all sequences
    std::cout << "nb_all_sequences after cms.generateMatrixCmsSeqs: gst.getNbAllSequences() " << gst.getNbAllSequences() << std::endl;
    
    return cms;
}

void save_common_motifs_matrix_to_csv(CommonMotifs cms, string output_csv_file){
    // according to test results, we will use saveMatrixCMS_ToCsv_File_dircrly
    /// cms.saveMatrixCMS_ToCsv_File_dircrly(output_csv_file);
    
    // print number of seqs, for debug:
    // we can't access cms->gst.getNbAllSequences(), because gst is private refernce in cms.
    // but when we access in cms.saveMatrixCMS_ToCsv_File, we get wrong answer 1.
    // when we access the same in cms.generateMatrixCmsSeqs() we get a correct answer.
    // it is the same object cms. I don't know where is the pbm.

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

}

int main(int argc, char *argv[]) {

    Args args = get_args(argc, argv);

    print_args(args);    

    // the output matrix in the csv file 

    string output_csv_file = generate_output_csv_file_name(args);
    cout << "C++ output_csv_file: " << output_csv_file << endl;
    
    auto start = chrono::high_resolution_clock::now();

    // sleep for 10 seconds, for debug only, to be removed
    sleep(3);

    //FastaFilesReader::getSaveInfosRNAFamiliesCSVFile(dir_name);

    //FastaFilesReader::groupALLFamiliesSeqsInOneFile(dir_name);

    // generate the common motifs matrix: desactivated for debug
    CommonMotifs cms = generate_common_motifs_matrix(args);
    
    auto end_before_saving = chrono::high_resolution_clock::now();

    cms.print_infos(); // for debug, need to reactivate this line of code

    // save the common motifs matrix to csv file
    save_common_motifs_matrix_to_csv(cms, output_csv_file); // desactivated for debug

    // only for debug:
    // generate random matrix of nb_families x 20 and save it in csv file output_csv_file
    /// generate_save_random_matrix(output_csv_file, args.nb_families);



    auto end = chrono::high_resolution_clock::now();

   
    // Calculating total time taken by the program.

    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end_before_saving - start).count();
    time_taken *= 1e-9;
    cout << "Time taken before save_to_csv is : " << fixed << time_taken << setprecision(9) << " sec" << endl;

    time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Time taken by whole program is : " << fixed << time_taken << setprecision(9) << " sec" << endl;


    return 0;
}