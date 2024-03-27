#include <iostream>
#include <string>
#include "Stockholm_SeedAlignment_To_plainFastaFiles.h"


// Struct to store the parsed command-line arguments
struct Args {
    std::string input_file_path; // -in
    std::string output_dir_path; // -out
};

void print_args_definition() {
    std::cout << "-in : <string> input file path for Stockholm file seed sequences." << std::endl;
    std::cout << "-out : <string> the output directory path" << std::endl;
}


void print_args(Args arg) {
    std::cout << "input file path: " << arg.input_file_path << std::endl;
    std::cout << "output directory path: " << arg.output_dir_path << std::endl;
}

// Function to parse command-line arguments using getopt
Args get_args(int argc, char* argv[]) {

    Args res;
 
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-in") {
            res.input_file_path = argv[++i];
        } else if (arg == "-out") {
            res.output_dir_path = argv[++i];
        } else {
            std::cout << "Usage: " << argv[0] << " -in <string> -out <string>" << std::endl;
            print_args_definition();
            exit(1);
        }
    }

    return res;
}


int main(int argc, char *argv[]) {

    Args args = get_args(argc, argv);

    std::cout << "file:  " << args.input_file_path << std::endl;

    uint32_t nb = Stockholm_SeedAlignment_To_plainFastaFiles::convert(args.input_file_path, args.output_dir_path);

    std::cout<<"Nb Total Families = "<<nb<<std::endl;

    return 0;
}

