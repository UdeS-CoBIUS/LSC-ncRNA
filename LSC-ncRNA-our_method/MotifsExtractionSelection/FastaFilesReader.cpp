//
// Created by ibra on 12/26/2019.
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include "FastaFilesReader.h"
#include "Utile.h"

#include <dirent.h>
#include <sys/stat.h>
#include <random>
#include <chrono>
#include <tuple>
#include <cstring>

#include <vector>
#include <string>



#ifdef _WIN32
#include <direct.h>
#endif

// must verify the end line to be \n
const char FastaFilesReader::newline = '\n'; // be careful on that, depending on the file \r\n, \r, \n ; you can see: https://stackoverflow.com/questions/6864759/determining-the-newline-character-for-the-environment-a-c-program-is-being-com

/*
vector<string> FastaFilesReader::getListSequences(const string &file_name)
{
    std::cout << "getListSequences --------- "<< std::endl;
    std::cout << "file_name: " << file_name << std::endl;

    vector<string> list_sequences;

    std::ifstream input(file_name);

    if(!input.good()){
        std::cout << " !input.good() -> return an empty vector" << std::endl;
        return list_sequences; // return an empty vector
    }

    std::cout << "input.good() " << std::endl;

    std::string line, name, content;

    while( std::getline( input, line ).good() )
    {
        if(line[0] == '>' ) // Identifier marker for a new entry
        {
            list_sequences.emplace_back("");
        }
        else
        {
            // maybe i can change this, to check the diffirent types of endline
            line.erase(std::remove(line.begin(),line.end(),newline),line.end());

            list_sequences.back() += line;
        }
    }

    input.close();

    std::cout << "list_sequences size: " << list_sequences.size() << std::endl;

    return list_sequences;
}
*/

vector<string> FastaFilesReader::getListSequences(const string &file_name)
{
    // get all the sequences in the file
    vector<string> list_sequences;

    std::ifstream input(file_name, std::ios::binary);  // Open in binary mode

    if (!input.good()) {
        std::cout << " !input.good() -> return an empty vector" << std::endl;
        return list_sequences; // return an empty vector
    }

    std::string line, sequence;
    char c;
    bool new_sequence = true;

    while (input.get(c)) {
        if (c == '>') {
            if (!sequence.empty()) {
                list_sequences.push_back(sequence);
                sequence.clear();
            }
            // Skip the rest of the header line
            std::getline(input, line);
            new_sequence = true;
        } else if (c != '\r' && c != '\n') {
            if (new_sequence) {
                sequence = c;
                new_sequence = false;
            } else {
                sequence += c;
            }
        }
    }

    // Add the last sequence if it exists
    if (!sequence.empty()) {
        list_sequences.push_back(sequence);
    }

    input.close();

    return list_sequences;
}



vector<vector<string>> FastaFilesReader::getListFamiliesSequences(const vector<string> &list_families_files_names)
{
    vector<vector<string>> list_families_sequences;
    list_families_sequences.reserve(list_families_files_names.size());

    for (const auto & file_name : list_families_files_names) {
        list_families_sequences.push_back(getListSequences(file_name));
    }

    return list_families_sequences;
}

inline char FastaFilesReader::separator()
{
#ifdef _WIN32
    return '\\';
#else
    return '/';
#endif
}


vector<vector<string>> FastaFilesReader::getListFamiliesSequences(const string &path_dir_families_files_names)
{
    // get all the files in the directory
    auto list_families_files_names = getListFiles(path_dir_families_files_names);
    
    auto list_families_sequences = getListFamiliesSequences(list_families_files_names);

    return list_families_sequences;
}


/*
 * Get File Name from a Path with or without extension
 * https://thispointer.com/c-how-to-get-filename-from-a-path-with-or-without-extension-boost-c17-filesytem-library/
 */
std::string FastaFilesReader::getFileName(const std::string filePath, bool withExtension = true, char seperator = FastaFilesReader::separator())
{
    // Get last dot position
    std::size_t dotPos = filePath.rfind('.');
    std::size_t sepPos = filePath.rfind(seperator);

    if(sepPos != std::string::npos)
    {
        return filePath.substr(sepPos + 1, filePath.size() - (withExtension || dotPos != std::string::npos ? 1 : dotPos) );
    }
    return "";
}


// I use C++14, for that i don't use the solution of (std::filesystem) of C++17
// https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
vector<string> FastaFilesReader::getListFiles(const string &path_dir)
{
    vector<string> list_files;

    DIR *dir;
    struct dirent *entry;
    struct stat file_stat{};

    dir = opendir (path_dir.c_str());

    if(dir== nullptr)
    {
        perror ("Could not open directory");
        return list_files; // return an empty list
    }

    string path_tmp;

    // get all the files within directory
    while ((entry = readdir (dir)) != nullptr)
    {
        // Skip hidden files and directories
        if (entry->d_name[0] == '.')
            continue;
        path_tmp=path_dir+FastaFilesReader::separator()+entry->d_name;

        if( stat(path_tmp.c_str(),&file_stat) <0) // get the stat
        {
            perror("stat error");
            continue; // go to the next entry
        }

        if (S_ISREG(file_stat.st_mode)) {
            //list_files.emplace_back(entry->d_name);
            list_files.emplace_back(path_tmp); // all the path not just the name of the file
        }
    }

    closedir (dir);

    std::sort(list_files.begin(),list_files.end());

    return list_files;
}


/*
// use c++17 // this looke like cause error...
vector<string> FastaFilesReader::getListFiles(const string &path_dir, const string &extension)
{
    vector<string> list_files;

    for (const auto & entry : fs::directory_iterator(path_dir))
    {
        if (fs::is_regular_file(entry) && entry.path().extension() == extension)
        {
            list_files.push_back(entry.path().string());
        }
    }

    std::sort(list_files.begin(), list_files.end());

    return list_files;
}
*/


vector<string> FastaFilesReader::get_Random_N_Files(const string &path_dir, unsigned int nb_files)
{
    auto list_families_files_names = FastaFilesReader::getListFiles(path_dir);

    // get a time-based seed
    unsigned seed = std::chrono::system_clock::now()
            .time_since_epoch()
            .count();

    auto rng = std::default_random_engine{seed};

    std::shuffle(list_families_files_names.begin(),list_families_files_names.end(),rng);

    vector<string>::const_iterator first = list_families_files_names.begin();
    vector<string>::const_iterator last = list_families_files_names.begin() + nb_files;

    vector<string> newVec(first, last);

    return newVec;
}

vector<string> FastaFilesReader::get_First_N_Files(const string &path_dir, unsigned int nb_files)
{
    auto list_families_files_names = FastaFilesReader::getListFiles(path_dir);

    vector<string>::const_iterator first = list_families_files_names.begin();
    vector<string>::const_iterator last = list_families_files_names.begin() + nb_files;

    vector<string> newVec(first, last);

    return newVec;
}

vector<vector<string>> FastaFilesReader::getListFamiliesSequences_FirstNFiles(const string &path_dir, uint32_t nb_files,
                                                                              uint32_t nb_seqs_allowed)
{
    auto list_families_files_names = FastaFilesReader::getListFiles(path_dir);

    uint32_t nb_families_added=0;

    vector<vector<string>> list_families_sequences;
    list_families_sequences.reserve(list_families_files_names.size());

    for (const auto & file_name : list_families_files_names)
    {
        if(isContainLessThanNSeqs(file_name,nb_seqs_allowed))
        {
            list_families_sequences.push_back(getListSequences(file_name));

            nb_families_added++;

            if(nb_families_added>nb_files) break;
        }
    }

    return list_families_sequences;
}


vector<vector<string>>
FastaFilesReader::getListFamiliesSequences_FirstNFiles_MinMax(const string &path_dir, uint32_t nb_files,
                                                              uint32_t nb_seqs_min, uint32_t nb_seqs_max)
{
    auto list_families_files_names = FastaFilesReader::getListFiles(path_dir);

    uint32_t nb_families_added=0;

    vector<vector<string>> list_families_sequences;
    list_families_sequences.reserve(nb_files);

    for (const auto & file_name : list_families_files_names)
    {
        if(FastaFilesReader::isNbSeqsBetween_MinMax(file_name, nb_seqs_min,nb_seqs_max))
        {
            list_families_sequences.push_back(getListSequences(file_name));

            nb_families_added++;

            if(nb_families_added>=nb_files) break;
        }
    }

    return list_families_sequences;
}

static void
get_Save_N_Random_Family_nbSeqs_in_MinMax(const string &path_dir_in, const string &path_dir_out, uint32_t nb_files,
                                                uint32_t nb_seqs_min, uint32_t nb_seqs_max)
{
    // get all files paths
    auto list_families_files_names = FastaFilesReader::getListFiles(path_dir_in);

    // get only files that have sequences between min and max
    vector<std::string> list_families_between_min_max;

    for (const auto & file_name : list_families_files_names)
    {
        if(FastaFilesReader::isNbSeqsBetween_MinMax(file_name, nb_seqs_min,nb_seqs_max))
        {
            list_families_between_min_max.push_back(file_name);
        }
    }

    // shuffle the list of families
    std::mt19937 rng(std::random_device{}());
    std::shuffle(list_families_between_min_max.begin(), list_families_between_min_max.end(), rng);

    // get the first n families
    vector<std::string> list_families_selected(list_families_between_min_max.begin(), list_families_between_min_max.begin() + nb_files);

    // save the files for the same name but to path_dir_out.
    for (const auto & file_name_path : list_families_selected)
    {
        std::string file_name_without_ext = FastaFilesReader::getFileName(file_name_path,false); 
        string path_writ = path_dir_out;
        path_writ += util::kPathSeparator; 
        path_writ += file_name_without_ext;
        path_writ += ".fasta";
        std::ifstream src(file_name_path, std::ios::binary);
        std::ofstream dst(path_writ, std::ios::binary);
        std::copy(std::istreambuf_iterator<char>(src),
                  std::istreambuf_iterator<char>(),
                  std::ostreambuf_iterator<char>(dst));

        // can use: 
        //  std::filesystem::copy_file(file_name_path, path_writ, std::filesystem::copy_options::overwrite_existing); // C++17
    }

    //-----------------------------------------
}

void FastaFilesReader::creat_dir_c(char* path_dir)
{
    struct stat st = {0};

    if (stat(path_dir, &st) == -1) {
        #ifdef _WIN32
            mkdir(path_dir);
        #else // Unix-like systems (Linux, macOS)
            mkdir(path_dir, 0777);
        #endif
    }
}

void
FastaFilesReader::construct_Train_Test_files(const string &path_dir_in, const string &path_dir_out, uint32_t nb_files,
                                             uint32_t nb_seqs_min, uint32_t nb_seqs_max,
                                             uint32_t percentage_nb_seqs_train)
{
    auto list_families_files_paths = FastaFilesReader::getListFiles(path_dir_in);

    uint32_t nb_families_added=0;

    uint32_t nb_seqs_train;
    uint32_t nb_seqs_file;
    //uint32_t nb_seqs_train = (nb_files*percentage_nb_seqs_train)/100;

    string buffer_train;
    string buffer_test;

    uint32_t idx;

    string train_dir = path_dir_out;
           train_dir += util::kPathSeparator;
           train_dir += "Train";

    string test_dir = path_dir_out;
           test_dir += util::kPathSeparator;
           test_dir += "Test";

    FastaFilesReader::creat_dir_c(const_cast<char *>(train_dir.c_str()));
    FastaFilesReader::creat_dir_c(const_cast<char *>(test_dir.c_str()));

    string file_name;

    for (const auto & file_path : list_families_files_paths)
    {
        nb_seqs_file = FastaFilesReader::getNbSeqs(file_path);
        if( nb_seqs_min<=nb_seqs_file  && nb_seqs_file <=nb_seqs_max)
        {
            nb_seqs_train = (nb_seqs_file*percentage_nb_seqs_train)/100;

            auto list_seqs=getListSequences(file_path);
            file_name = FastaFilesReader::getFileName(file_path);

            for (idx = 0; idx < nb_seqs_train; ++idx)
            {
                buffer_train += ">";
                buffer_train += file_name;
                buffer_train += "_";
                buffer_train += util::to_string(idx);
                buffer_train += "\n";
                buffer_train += list_seqs.at(idx);
                buffer_train += "\n";
            }

            write_to_file(buffer_train, file_name, train_dir);
            buffer_train.clear();

            for (; idx < nb_seqs_file; ++idx)
            {
                buffer_test += ">";
                buffer_test += file_name;
                buffer_test += "_";
                buffer_test += util::to_string(idx);
                buffer_test += "\n";
                buffer_test += list_seqs.at(idx);
                buffer_test += "\n";
            }

            write_to_file(buffer_test, file_name, test_dir);
            buffer_test.clear();

            nb_families_added++;

            if(nb_families_added>=nb_files) break;
        }
    }
}

void FastaFilesReader::write_to_file(const string &buffer, const string &file_name, const string &dir)
{

    ofstream file_output;

    string path_writ = dir;
           path_writ += util::kPathSeparator;
           path_writ += file_name;

    cout<<"path : "<<path_writ<<endl;

    file_output.open(path_writ);

    file_output<<buffer;

    file_output.close();

}

bool FastaFilesReader::isContainLessThanNSeqs(const string &file_name, uint32_t nb_seqs_allowed)
{

    std::ifstream input(file_name);

    if(!input.good())
        return false; // return an empty vector

    std::string line, name, content;
    uint32_t nb=0;

    while( std::getline( input, line ).good() )
    {
        if(line[0] == '>' ) // Identifier marker for a new entry
        {
            nb++;

            if(nb>nb_seqs_allowed) return false;
        }
    }

    input.close();

    return true;
}

/**
 *
 * @param file_name
 * @param min
 * @param max
 * @return a list of files names, where each file contain min<=nb-of-seq<=max
 */
bool FastaFilesReader::isNbSeqsBetween_MinMax(const string &file_name, uint32_t min, uint32_t max)
{
    std::ifstream input(file_name);

    if(!input.good())
        return false; // return an empty vector

    std::string line, name, content;
    uint32_t nb=0;

    while( std::getline( input, line ).good() )
    {
        if(line[0] == '>' ) // Identifier marker for a new entry
        {
            nb++;

            if(nb>max) return false;
        }
    }

    input.close();

    return nb>=min;
}


int FastaFilesReader::getNbSeqs(const string &file_name)
{
    std::ifstream input(file_name);

    if(!input.good())
        return false; // return an empty vector

    std::string line, name, content;
    uint32_t nb=0;

    while( std::getline( input, line ).good() )
    {
        if(line[0] == '>' ) // Identifier marker for a new entry
        {
            nb++;
        }
    }

    input.close();

    return nb;
}


string FastaFilesReader::getNbMinMaxMeanFamilySequences(const string &file_name)
{
    uint32_t nb_seqs=0;
    uint32_t min_seq = 0xffffffff; // max value : 4'294'967'295
    uint32_t max_seq=0;
    uint32_t all_size = 0;
    string family_info;

    std::ifstream input(file_name);

    if(!input.good())
        return "could not open the file("+file_name+")...."; // return : nb=0, min=0,max0

    std::string line, sequence;

    bool is_first_sequence_readed = false;

    while( std::getline( input, line ).good() )
    {
        if(line[0] == '>' ) // Identifier marker for a new entry
        {
            nb_seqs++;

            if(!is_first_sequence_readed)
            {
                is_first_sequence_readed = true;

            } else
            {
                if(sequence.size() > max_seq) max_seq = sequence.size();

                if(sequence.size() < min_seq) min_seq = sequence.size();

                all_size += sequence.size();
                sequence.clear();
            }

        }
        else
        {
            // maybe i can change this, to check the diffirent types of endline
            line.erase(std::remove(line.begin(),line.end(),newline),line.end());

            sequence += line;
        }
    }

    // treat the last sequence
    if(sequence.size() > max_seq) max_seq = sequence.size();
    if(sequence.size() < min_seq) min_seq = sequence.size();
    all_size += sequence.size();

    input.close();

    family_info+=file_name;
    family_info+=CSV_SEPARATOR;
    family_info+=util::to_string(nb_seqs);
    family_info+=CSV_SEPARATOR;
    family_info+=util::to_string(min_seq);
    family_info+=CSV_SEPARATOR;
    family_info+=util::to_string(max_seq);
    family_info+=CSV_SEPARATOR;

    if (nb_seqs > 0) {
        family_info += util::to_string(all_size / nb_seqs);
    } else {
        family_info += "0";
    } // mean
    family_info+=CSV_SEPARATOR;

    return family_info;
}

vector<string> FastaFilesReader::getInfosListFamiliesSequences(const string &path_dir)
{
    auto list_families_files_names = FastaFilesReader::getListFiles(path_dir);

    vector<string> list_families_infos;
    list_families_infos.reserve(list_families_files_names.size());


    for (const auto & family_file_name : list_families_files_names)
    {
        list_families_infos.push_back(getNbMinMaxMeanFamilySequences(family_file_name));
    }

    return list_families_infos;
}

void FastaFilesReader::getSaveInfosRNAFamiliesCSVFile(const string &path_dir)
{
    auto list_infos_RNAFamilies = getInfosListFamiliesSequences(path_dir);

    string buffer;
    ofstream file_csv;

    for(const auto &family_info:list_infos_RNAFamilies)
    {
        buffer+=family_info;
        buffer+="\n";
    }

    file_csv.open ("All_RNA_Families_infos.csv");
    cout<<" file_csv.open (\"All_RNA_Families_infos.csv\"); "<<endl;

    file_csv<<buffer;

    file_csv.close();
    cout<<"file_csv.close();"<<endl;
}

void FastaFilesReader::groupALLFamiliesSeqsInOneFile(const string &path_dir)
{
    auto list_families_sequences = FastaFilesReader::getListFamiliesSequences(path_dir);

    string buffer;
        buffer+="id_global_seq";
        buffer+=CSV_SEPARATOR;
        buffer+="id_family";
        buffer+=CSV_SEPARATOR;
        buffer+="sequence";
        buffer+=newline;

    uint32_t id_seq_total=0; // idx of seqs in all database
    uint32_t id_family=0; // the family id, RFxxxxx , we use only xxxxx

    for(const auto &list_family:list_families_sequences)
    {
        for (const auto & seq : list_family)
        {
            buffer+=util::to_string(id_seq_total);
            buffer+=CSV_SEPARATOR;
            buffer+=util::to_string(id_family);
            buffer+=CSV_SEPARATOR;
            buffer+=seq;
            buffer+=newline;

            id_seq_total++;
        }

        id_family++;
    }

    cout << "total : " << id_seq_total << endl;


    ofstream file_csv;

    file_csv.open ("All_RNA_Families_IdTotal_IdFamily_Seq.csv");
    cout<<" file_csv.open (\"All_RNA_Families_IdTotal_IdFamily_Seq.csv\"); "<<endl;

    file_csv<<buffer;

    file_csv.close();
    cout<<"file_csv.close();"<<endl;

}
