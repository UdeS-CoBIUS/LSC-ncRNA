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
    // by default:
    //string dir_input = R"(C:\Users\ibra\OneDrive - USherbrooke\DatatSet\Rfam_RNA_Seed_ungapedFasta)"; // this not used.
	// this next the based files used as Rfam 14.1 in our tests.
    //string dir_input = R"(C:\Users\ibra\OneDrive - USherbrooke\DatatSet\Rfam_RNAFamilies_Stockholm_SeedAlignment_To_PlainFastaFiles\RNA_Family)";
    string dir_input = R"(C:\Users\ibra\Desktop\Infernal\Secondary_and_not_Rfam_some_family_v14.7\SecondaryBased_lessSimilarSeqs)";
    string dir_output = R"(C:\Users\ibra\Desktop\Infernal\Secondary_and_not_Rfam_some_family_v14.7\SecondaryBased_lessSimilarSeqs_train_test)";

    uint32_t nb_families = 600;
    uint32_t min_nb_seqs_allowed = 20;
    uint32_t max_nb_seqs_allowed = 105;

    uint32_t percentage_nb_seqs_train = 70;


    //getMainArgv(argc, argv, dir_input, dir_output, nb_families, min_nb_seqs_allowed, max_nb_seqs_allowed, percentage_nb_seqs_train);

    cout << "dir input : " << dir_input << endl;
    cout << "dir output : " << dir_output << endl;

    cout << "nb families : " << nb_families << endl;
    cout << "Min nb seqs : " << min_nb_seqs_allowed << endl;
    cout << "Max nb seqs : " << max_nb_seqs_allowed << endl;
    cout << "percentage_nb_seqs_train : " << percentage_nb_seqs_train << endl;


	//FastaFilesReader::getSaveInfosRNAFamiliesCSVFile(dir_input); // this used to get all informatio as nb seq, min seq len, max seq len, average seq len , and save to csv file.


//    FastaFilesReader::construct_Train_Test_files(dir_input, dir_output, nb_families,
//                                                min_nb_seqs_allowed, max_nb_seqs_allowed,
//                                                percentage_nb_seqs_train);

    FastaFilesReader::construct_Train_Test_files(dir_input, dir_output,
                                                 min_nb_seqs_allowed,
                                                 percentage_nb_seqs_train);

    //FastaFilesReader::get_Families_files(dir_input, dir_output, 1000,
    //                                            min_nb_seqs_allowed, max_nb_seqs_allowed);



//    FastaFilesReader::copy_Train_Test_files(dir_input, dir_output, nb_families,
//                                                 min_nb_seqs_allowed, max_nb_seqs_allowed);

    return 0;
}