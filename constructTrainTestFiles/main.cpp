#include <iostream>
#include "FastaFilesReader.h"

void getMainArgv(int argc, char *argv[],
                string &dir_input, string &dir_output, uint32_t &nb_families,
                 uint32_t &min_nb_seqs_allowed,uint32_t & max_nb_seqs_allowed,
                 uint32_t &percentage_nb_seqs_train)
{

    uint32_t nb_arg = 6;

    if(argc==nb_arg+1)
    {
        dir_input = argv[1];
        dir_output = argv[2];

        nb_families = strtol(argv[3], nullptr, 10);
        min_nb_seqs_allowed = strtol(argv[4], nullptr, 10);
        max_nb_seqs_allowed = strtol(argv[5], nullptr, 10);

        percentage_nb_seqs_train = strtol(argv[6], nullptr, 10);
    }


    /*
    // use getopt
    int c;
    extern char *optarg;
    extern int optind, optopt;

    while ((c = getopt(argc, argv, ":ionmxp")) != -1) {
        switch(c) {
            case 'i': dir_input = optarg; break;
            case 'o': dir_output = optarg; break;
            case 'n': nb_families = strtol(optarg, nullptr, 10); break;
            case 'm': min_nb_seqs_allowed = strtol(optarg, nullptr, 10); break;
            case 'x': max_nb_seqs_allowed = strtol(optarg, nullptr, 10); break;
            case 'p': percentage_nb_seqs_train = strtol(optarg, nullptr, 10); break;
            case ':':       /// -f or -o without operand
                fprintf(stderr, "Option -%c requires an operand\n", optopt);
                break;
            case '?':
                fprintf(stderr,
                        "Unrecognized option: -%c\n", optopt);
            default: fprintf(stderr,
                             "Error: -%c\n", optopt);
        }
    }

    */
}

int main(int argc, char *argv[])
{
    string dir_input = R"(C:\Users\ibra\Desktop\Infernal\Rfam_14.1_dataset\Rfam_more_than_3seqs\Train)";
    string dir_output = R"(C:\Users\ibra\Desktop\Infernal\Rfam_14.1_dataset\Rfam_350)";

    uint32_t nb_families = 165; // just to use all families
    uint32_t min_nb_seqs_allowed = 4;
    uint32_t max_nb_seqs_allowed = 10; // max is 1542 i think

     dir_input = argv[1];
     dir_output = argv[2];
     nb_families = std::stoi(argv[3]);
     min_nb_seqs_allowed = std::stoi(argv[4]);
     max_nb_seqs_allowed = std::stoi(argv[5]);

    /// uint32_t percentage_nb_seqs_train = 70;// change to use percentage_nb_seqs_test
    uint32_t percentage_nb_seqs_test = 30; // 100 - percentage_nb_seqs_train
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


    //getMainArgv(argc, argv, dir_input, dir_output, nb_families, min_nb_seqs_allowed, max_nb_seqs_allowed, percentage_nb_seqs_train);

    cout << "dir input : " << dir_input << endl;
    cout << "dir output : " << dir_output << endl;

    cout << "nb families : " << nb_families << endl;
    cout << "Min nb seqs : " << min_nb_seqs_allowed << endl;
    cout << "Max nb seqs : " << max_nb_seqs_allowed << endl;
    cout << "percentage_nb_seqs_train : " << 100 - percentage_nb_seqs_test << endl;
    cout << "percentage_nb_seqs_test : " << percentage_nb_seqs_test << endl;



	//FastaFilesReader::getSaveInfosRNAFamiliesCSVFile(dir_input); // this used to get all informatio as nb seq, min seq len, max seq len, average seq len , and save to csv file.


    // FastaFilesReader::construct_Train_Test_files(dir_input, dir_output, nb_families,
    //                                            min_nb_seqs_allowed, max_nb_seqs_allowed,
    //                                            percentage_nb_seqs_test);

    FastaFilesReader::get_Save_N_Random_Family_nbSeqs_in_MinMax(dir_input,dir_output,nb_families,min_nb_seqs_allowed,max_nb_seqs_allowed);

//    FastaFilesReader::construct_Train_Test_files(dir_input, dir_output,
//                                                 min_nb_seqs_allowed,
//                                                 percentage_nb_seqs_test);

    //FastaFilesReader::get_Families_files(dir_input, dir_output, 1000,
    //                                            min_nb_seqs_allowed, max_nb_seqs_allowed);



//    FastaFilesReader::copy_Train_Test_files(dir_input, dir_output, nb_families,
//                                                 min_nb_seqs_allowed, max_nb_seqs_allowed);

    return 0;
}