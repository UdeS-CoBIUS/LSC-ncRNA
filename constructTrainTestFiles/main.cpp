#include <iostream>
#include "FastaFilesReader.h"

// Struct to store the parsed command-line arguments
struct Args {
    string dir_input; // -in
    string dir_output; // -out
    string mode = "-i"; // -m
    uint32_t nb_families; // -nf
    uint32_t min_nb_seqs_allowed = 4; // -mins
    uint32_t max_nb_seqs_allowed = 1000; // -maxs
    uint32_t percentage_nb_seqs_test = 30; // -pt

    /// uint32_t percentage_nb_seqs_train = 70;// change to use percentage_nb_seqs_test
    // uint32_t percentage_nb_seqs_test = 30; // 100 - percentage_nb_seqs_train
    // we use percentage_nb_seqs_test to compute the number of test sequences
    // the remaining are dirctly the train sequences.
    // And we want give the train more sequences that the test
    // because in our old code: we use (nb_seqs_train = (nb_seqs_file*percentage_nb_seqs_train)/100;)
    // se we use the devision function
    // 4*70/100 = 2  (2.8, but the devision give the integer value wich is 2)
    // so , number train sequences is 2, and the remaining is for test which is also 2.
    // and this is a problem, because like this we have 50% 50%
    // se begin by test sequences like
    // (nb_seqs_test = (nb_seqs_file*percentage_nb_seqs_test)/100;)
    // like this we will have less sequences in test part, because of the division.
    // 4*30/100 = 1 ( real result is 1.2), and remaing 3 are for train part.

};

void print_args_definition() {
    cout << "-in : <string> a path directory for fasta files" << endl;
    cout << "-out : <string> a main path directory the results out" << endl;
    cout << "-nf : <integer> number of families "<< endl;
    cout << "-mins : <integer>, min number of sequences (default 4)"<< endl;
    cout << "-maxs : <integer>, max number of sequences"<< endl;
    cout << "-pt : <integer>, percentage number sequences Test" << endl;

    cout << "-m : <string> the used mode, several mode are avaialable: "<< endl;
    cout << "-m i: informatio, get all informations as nb seqs, min seq len, max seq len, average seq len , and save to csv file." << endl;
    cout << "-m s: Sample dataset, get n random families that have nb seqs between min and max, and save them to dir_output" << endl;
    cout << "-m sttmm : Split Tarin Test Min Max, for a given nb of families, and min max number of seqs, split to train and test" << endl;
    cout << "-m sttm : Split Tarin Test Min, consider only min number of seqs, and split all files in input folder to train and test" << endl;
    cout << "-m stt : Split Tarin Test, for all files in input folder split to train and test" << endl;
}


void print_args(Args arg) {
    cout << "dir_input: " << arg.dir_input << endl;
    cout << "dir_output: " << arg.dir_output << endl;
    cout << "nb_families: " << arg.nb_families << endl;
    cout << "min_nb_seqs_allowed: " << arg.min_nb_seqs_allowed << endl;
    cout << "max_nb_seqs_allowed: " << arg.max_nb_seqs_allowed << endl;
    cout << "percentage_nb_seqs_test: " << arg.percentage_nb_seqs_test << endl;
    cout << "mode: " << arg.mode << endl;
}

// Function to parse command-line arguments
Args get_args(int argc, char* argv[]) {

    Args res;
 
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-in") {
            res.dir_input = argv[++i];
        } else if (arg == "-out") {
            res.dir_output = argv[++i];
        } else if (arg == "-nf") {
            res.nb_families = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-mins") {
            res.min_nb_seqs_allowed = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-maxs") {
            res.max_nb_seqs_allowed = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-pt") {
            res.percentage_nb_seqs_test = strtol(argv[++i], nullptr, 10);
        } else if (arg == "-m") {
            res.mode =argv[++i];
        } else {
            cout << "Usage: " << argv[0] << " -in <string> -out <string> -nf <integer> -mins <integer> -maxs <integer> -pt <integer> -m <string>" << endl;
            print_args_definition();
            exit(1);
        }
    }

    return res;
}

int main(int argc, char *argv[])
{

    Args args = get_args(argc, argv);
    
    print_args(args);
    

    if(args.mode == "i"){ // -m i: information
        // this used to get all informatio as nb seq, min seq len, max seq len, average seq len , and save to csv file.
        cout << " -m i: get and Save to csv all ncRNAs families informations" << endl;
	    FastaFilesReader::getSaveInfosRNAFamiliesCSVFile(dir_input); 
    }
    else if(args.mode=="s"){ // -m s: sample
        // this for get N random families that have nb seqs between min and max, and save them to dir_output
        cout << " -m s: get N random families that have nb seqs between min and max, and save them to dir_output..." << endl;
        FastaFilesReader::get_Save_N_Random_Family_nbSeqs_in_MinMax(dir_input,dir_output,nb_families,min_nb_seqs_allowed,max_nb_seqs_allowed);
    }
    else if(args.mode == "sttmm"){ // -m sttmm : split tarin test min max
        cout << " -m sttmm : split tarin test min max ..." << endl;
        FastaFilesReader::construct_Train_Test_files(dir_input, dir_output, nb_families, min_nb_seqs_allowed, max_nb_seqs_allowed, percentage_nb_seqs_test);
    }
    else if(args.mode == "sttm"){ // -m sttm : split tarin test min
        cout << " -m sttm : split tarin test min ..." << endl;
        FastaFilesReader::construct_Train_Test_files(dir_input, dir_output, min_nb_seqs_allowed, percentage_nb_seqs_test);
    }
    else if(args.mode == "stt"){ // -m stt : split tarin test
        // split the all files in the directory into Train Test files.
        cout << " -m stt : split tarin test ..." << endl;
        FastaFilesReader::construct_Train_Test_files(dir_input, dir_output, percentage_nb_seqs_test);
    }
    
    //FastaFilesReader::get_Families_files(dir_input, dir_output, 1000,
    //                                            min_nb_seqs_allowed, max_nb_seqs_allowed);



//    FastaFilesReader::copy_Train_Test_files(dir_input, dir_output, nb_families,
//                                                 min_nb_seqs_allowed, max_nb_seqs_allowed);

    return 0;
}