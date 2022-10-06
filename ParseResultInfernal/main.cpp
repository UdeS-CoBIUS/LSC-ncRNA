#include <iostream>
#include <fstream>
#include <string>

#include <filesystem>
#include <iterator>
#include "FastaFilesReader.h"
#include "Utile.h"
#include <cmath>

namespace fs = std::filesystem;

/**
 *  split string with white space
 *  form : https://stackoverflow.com/questions/2275135/splitting-a-string-by-whitespace-in-c
 * @param input
 * @return
 */
std::vector<std::string> split(std::string const &input)
{
    std::istringstream buffer(input);

    std::vector<std::string> ret{std::istream_iterator<std::string>(buffer),
                                 std::istream_iterator<std::string>()};
    return ret;
}


std::string getFileName_without_extenstion(const std::string file_with_extension)
{
    // Get last dot position
    std::size_t dotPos = file_with_extension.find('.');


    return file_with_extension.substr(0,dotPos);
}

std::vector<std::string> classParser(const std::string &file_name)
{
    std::vector<std::string> list_infos;

    std::ifstream input(file_name);

    if(!input.good())
        return list_infos; // return an empty vector

    std::string query_marker ="Query:";
    std::string res_found__marker ="  (";
    std::string res_found_1__marker ="  (1)";
    std::string no_hits_detected__marker ="   [No hits detected that satisfy reporting thresholds]";
    std::string only_hmm_marker ="Total HMM hits reported:";
    std::string cm_marker ="Total CM hits reported:"; // it include the hmm filter in its pipline
    std::string total_cm_and_hmm_marker ="Total CM and HMM hits reported:";

    std::string line;
    std::string query_entry;
    std::string query_entry_family;
    std::string content;

    uint32_t nb_entry_total=0;
    uint32_t nb_entry_no_hits_detected=0;

    uint32_t nb_entry_correct_class=0;
    uint32_t nb_entry_incorrect_class=0;
    uint32_t nb_entry_other_class_from_2_res=0;

    uint32_t nb_only_hmm_hits_found = 0;
    uint32_t nb_cm_hits_found = 0;
    uint32_t nb_total_cm_and_hmm_hits_found = 0;

    while( std::getline( input, line ).good() )
    {
        if (line.rfind(query_marker, 0) == 0) // Identifier marker for a new entry
        {
            nb_entry_total++;
            query_entry = split(line).at(1);

            query_entry_family = getFileName_without_extenstion(query_entry);

            list_infos.emplace_back( query_entry );
        }
        else
        if (line.rfind(res_found_1__marker, 0) == 0
           &&
           split(line).at(5)!="cm" && split(line).at(5)!="hmm" ) // the firs result
        {
            list_infos.back() += line;

            if(split(line).at(5).find(query_entry_family) != std::string::npos)
            {
                nb_entry_correct_class++;
            }
            else
            {
                nb_entry_incorrect_class++;

                std::cout << "query_entry_family = "<<query_entry_family<< " | line : " << line << std::endl;
            }
        }
        else
        if (line.rfind(res_found__marker, 0) == 0
        &&
        split(line).at(5)!="cm" && split(line).at(5)!="hmm" ) // other result from 2 and so on
        {
            list_infos.back() += line;

            nb_entry_other_class_from_2_res++;
        }
        else
        if (line.rfind(no_hits_detected__marker, 0) == 0 ) // Identifier marker for a new entry
        {
            nb_entry_no_hits_detected++;
            list_infos.back() += line;
        }
        else
        if (line.rfind(only_hmm_marker, 0) == 0)
        {
            //std::cout<<"only_hmm_marker : "<<line <<" | "<<split(line).at(4)<<" | "<<strtol(split(line).at(4).c_str(), nullptr, 10)<<std::endl;
            nb_only_hmm_hits_found += strtol(split(line).at(4).c_str(), nullptr, 10);
        }
        else
        if (line.rfind(cm_marker, 0) == 0)
        {
            nb_cm_hits_found += strtol(split(line).at(4).c_str(), nullptr, 10);
        }
        else
        if (line.rfind(total_cm_and_hmm_marker, 0) == 0)
        {
            nb_total_cm_and_hmm_hits_found += strtol(split(line).at(6).c_str(), nullptr, 10);
        }
    }

    input.close();

    std::cout << "nb_entry_total = " << nb_entry_total << std::endl;
    std::cout << "nb_entry_no_hits_detected  = " << nb_entry_no_hits_detected/2 << std::endl;
    std::cout << "nb_entry_correct_class  = " << nb_entry_correct_class << std::endl;
    std::cout << "nb_entry_incorrect_class  = " << nb_entry_incorrect_class << std::endl;
    std::cout << "nb_entry_other_class_from_2_res  = " << nb_entry_other_class_from_2_res << std::endl;

    std::cout << "nb_only_hmm_hits_found  = " << nb_only_hmm_hits_found << std::endl;
    std::cout << "nb_cm_hits_found  = " << nb_cm_hits_found << std::endl;
    std::cout << "nb_total_cm_and_hmm_hits_found  = " << nb_total_cm_and_hmm_hits_found << std::endl;

    return list_infos;
}


/**
 *
 * @param file_name
 * @return   list (family_name, list(begin,end))  from the file pars
 */
std::vector< std::pair<std::string, std::vector<std::pair<uint32_t , uint32_t >> > >  SearchParser(const std::string &file_name)
{
    //std::vector<std::string> list_infos;

    std::vector< std::pair<std::string, std::vector<std::pair<uint32_t , uint32_t >> > > list_infos;
    //string
    // vector: (begin end)

    std::ifstream input(file_name);

    if(!input.good())
        return list_infos; // return an empty vector

    std::string query_marker ="Query:";
    std::string res_found__marker ="  (";
    std::string no_hits_detected__marker ="   [No hits detected that satisfy reporting thresholds]";
    std::string only_hmm_marker ="Total HMM hits reported:";
    std::string cm_marker ="Total CM hits reported:"; // it include the hmm filter in its pipline
    std::string total_cm_and_hmm_marker ="Total CM and HMM hits reported:";

    std::string line;
    std::string query_entry;
    std::string query_entry_family;
    std::string content;

    uint32_t nb_entry_total=0;
    uint32_t nb_entry_no_hits_detected=0;

    uint32_t nb_entry_correct_class=0;
    uint32_t nb_entry_incorrect_class=0;
    uint32_t nb_entry_other_class_from_2_res=0;

    uint32_t nb_only_hmm_hits_found = 0;
    uint32_t nb_cm_hits_found = 0;
    uint32_t nb_total_cm_and_hmm_hits_found = 0;

    while( std::getline( input, line ).good() )
    {
        if (line.rfind(query_marker, 0) == 0) // Identifier marker for a new entry
        {
            nb_entry_total++;
            query_entry = split(line).at(1);

            query_entry_family = getFileName_without_extenstion(query_entry);

            query_entry_family = query_entry_family.substr(query_entry_family.find("RF"));

            list_infos.emplace_back( std::make_pair(query_entry_family, std::vector<std::pair<uint32_t , uint32_t >>()) );
        }
        else
        //if (line.rfind(res_found__marker, 0) == 0
        if ( !split(line).empty() && split(line).at(0).rfind('(',0)==0
            &&
            split(line).at(5)!="cm" && split(line).at(5)!="hmm" )
        {
            uint32_t start =  strtol(split(line).at(6).c_str(), nullptr, 10);
            uint32_t end   =  strtol(split(line).at(7).c_str(), nullptr, 10);

            list_infos.back().second.emplace_back(std::make_pair(start,end));
        }
        else
        if (line.rfind(no_hits_detected__marker, 0) == 0 ) // Identifier marker for a new entry
        {
            nb_entry_no_hits_detected++;
        }
        else
        if (line.rfind(only_hmm_marker, 0) == 0)
        {
            nb_only_hmm_hits_found += strtol(split(line).at(4).c_str(), nullptr, 10);
        }
        else
        if (line.rfind(cm_marker, 0) == 0)
        {
            nb_cm_hits_found += strtol(split(line).at(4).c_str(), nullptr, 10);
        }
        else
        if (line.rfind(total_cm_and_hmm_marker, 0) == 0)
        {
            nb_total_cm_and_hmm_hits_found += strtol(split(line).at(6).c_str(), nullptr, 10);
        }
    }

    input.close();

    std::cout << "nb_entry_total = " << nb_entry_total << std::endl;
    std::cout << "nb_entry_no_hits_detected  = " << nb_entry_no_hits_detected/2 << std::endl;
    //std::cout << "nb_entry_correct_class  = " << nb_entry_correct_class << std::endl;
    //std::cout << "nb_entry_incorrect_class  = " << nb_entry_incorrect_class << std::endl;
    //std::cout << "nb_entry_other_class_from_2_res  = " << nb_entry_other_class_from_2_res << std::endl;

    std::cout << "nb_only_hmm_hits_found  = " << nb_only_hmm_hits_found << std::endl;
    std::cout << "nb_cm_hits_found  = " << nb_cm_hits_found << std::endl;
    std::cout << "nb_total_cm_and_hmm_hits_found  = " << nb_total_cm_and_hmm_hits_found << std::endl;



    for (auto & k : list_infos) {
        std::sort(k.second.begin(), k.second.end());
    }

    return list_infos;
}

std::vector< std::pair< std::string, std::vector< std::pair<uint32_t , uint32_t > > > >
        getStartENd_In_FakeGenome_single_char_and_with_RealRNASeqs(size_t n, const std::vector< std::pair<std::string, std::vector<std::string> > > & list_RNASeq_byFamilies, char nuc='A')
{

    std::vector< std::pair< std::string, std::vector< std::pair<uint32_t , uint32_t > > > > list_seqs_start_end_by_families;
    //list_seqs_start_end_by_families.first = family name;
    //list_seqs_start_end_by_families.second = vecor of pair(start,end);

    uint32_t nb_total_seqs = 0;
    uint32_t size_all_RNA_together = 0;

    for (const auto & list_RNASeq : list_RNASeq_byFamilies)
    {
        for (const auto & rna : list_RNASeq.second) {
            size_all_RNA_together+=rna.size();
            nb_total_seqs++;
        }
    }

    cout<<" nb_total_seqs = "<<nb_total_seqs <<endl;
    cout<<" size_all_RNA_together = "<<size_all_RNA_together <<endl;


    //return list_seqs_start_end_by_families;


    uint32_t nb_char_for_random_genome;
    nb_char_for_random_genome = (n > size_all_RNA_together) ? (n - size_all_RNA_together) : 0;

    cout<<" nb_char_for_random_genome = "<<nb_char_for_random_genome <<endl;

    std::string random_genome(nb_char_for_random_genome,nuc);

    std::string genome_with_real_RNA;
    genome_with_real_RNA.reserve(n);

    uint32_t offset_insert = nb_char_for_random_genome/nb_total_seqs;

    cout<<" offset_insert = "<< offset_insert  <<endl;

    uint32_t start = 0;
    uint32_t end = 0;


    string family_name;
    std::vector< std::pair<uint32_t , uint32_t > > vect_start_end;


    uint32_t i = 0;
    for (const auto & list_RNASeq : list_RNASeq_byFamilies)
    {
        family_name = list_RNASeq.first;

        for (const auto & rna : list_RNASeq.second)
        {
            start = genome_with_real_RNA.length();
            genome_with_real_RNA += rna;
            end = genome_with_real_RNA.length();

            genome_with_real_RNA += random_genome.substr(offset_insert*i ,offset_insert);

            vect_start_end.emplace_back(std::make_pair(start,end));

            i++;
        }

        list_seqs_start_end_by_families.emplace_back(family_name,vect_start_end);

        vect_start_end.clear();
    }


    genome_with_real_RNA += random_genome.substr(offset_insert*i);


    for (auto & k : list_seqs_start_end_by_families) {
        std::sort(k.second.begin(), k.second.end());
    }

    return list_seqs_start_end_by_families;
}


std::vector< std::pair<std::string, std::vector<std::string> > >
        getList_RNASeqs_by_family(const std::string& dir_test)
{
    auto list_families_files_names = FastaFilesReader::getListFiles(dir_test);
    auto list_families_sequences = FastaFilesReader::getListFamiliesSequences(list_families_files_names);


    std::vector< std::pair<std::string, std::vector<std::string> > > list_RNASeqs_by_family;

    std::string family_name;

    for (uint32_t i = 0; i < list_families_files_names.size(); ++i)
    {
        family_name = getFileName_without_extenstion(FastaFilesReader::getFileName(list_families_files_names.at(i)));

        list_RNASeqs_by_family.emplace_back( std::make_pair(family_name, list_families_sequences.at(i)) );
    }

    return list_RNASeqs_by_family;
}
void check_search_result(const string &dir_test, const string &file_to_parse)
{

    // -----------------------------------------------------------------------------------------
    // 1 ) get start end poition of RNAs by families in the genrated fake genome
    // -----------------------------------------------------------------------------------------

    auto list_RNASeqs_by_family = getList_RNASeqs_by_family(dir_test);

    uint32_t n = 3'000'000;
    auto list_start_end_in_genome = getStartENd_In_FakeGenome_single_char_and_with_RealRNASeqs(n,list_RNASeqs_by_family,'A');

    int offest = list_start_end_in_genome.at(0).second.at(1).first
                      -
                        list_start_end_in_genome.at(0).second.at(0).second;

    cout<<" offest = " << offest <<endl;

    for (auto & k : list_start_end_in_genome) {
        cout<<k.first<<" : ----------------------------- \n";
        for (auto & xy : k.second)
        {
            cout << xy.first << " , "<<xy.second << endl;
        }
    }


    // -----------------------------------------------------------------------------------------
    // 2) parse the result file of cmsearch and get the resilts of start end of the found RNA by mamilies
    // -----------------------------------------------------------------------------------------
    auto list_res_parsedFile = SearchParser(file_to_parse);

    for (auto & k : list_res_parsedFile) {
        cout<<k.first<<" : ----------------------------- \n";
        for (auto & xy : k.second)
        {
            cout << xy.first << " , "<<xy.second << endl;
        }
    }

    // -----------------------------------------------------------------------------------------
    // 3 ) check the resulta from the parsed result file with the real genrated genome postoion of RNA by families
    // -----------------------------------------------------------------------------------------

    uint32_t inter_f_start = 0;
    uint32_t inter_f_end = 0;

    uint32_t nb_entry_correct_class=0;
    uint32_t nb_entry_incorrect_class=0;
    uint32_t nb_entry_other_class_from_2_res=0;

    for (uint32_t i = 0; i < list_res_parsedFile.size(); ++i)
    {
        if(list_res_parsedFile.at(i).first == list_start_end_in_genome.at(i).first)
        {
            cout << "fam = "<<list_res_parsedFile.at(i).first;

            inter_f_start = list_start_end_in_genome.at(i).second.at(0).first;
            inter_f_end = list_start_end_in_genome.at(i).second.back().second;

            cout<<"[ "<<inter_f_start << " , "<<inter_f_end <<"] "<<endl;

            for (auto & j : list_res_parsedFile.at(i).second)
            {
                //cout<<"("<<j.first<<","<<j.second<<") ; ";
                if(j.second < inter_f_start || j.first > inter_f_end )
                {
                    nb_entry_incorrect_class++;
                }
                else
                {
                    // NAIF TEST for all element in the families

                    const uint32_t found_seq_start = j.first;
                    const uint32_t found_seq_end = j.second;

                    uint32_t real_seq_start = 0;
                    uint32_t real_seq_end = 0;
                    uint32_t real_seq_length = 0;

                    bool found = false;

                    for (auto & k : list_start_end_in_genome.at(i).second)
                    {
                        real_seq_start = k.first;
                        real_seq_end = k.second;
                        real_seq_length = real_seq_end - real_seq_start;

                        //int x = (labs(real_seq_end - found_seq_end) - labs(real_seq_start - found_seq_start));
                        int x = labs(real_seq_end - found_seq_end);

                        if(x<real_seq_length) // condition wrong, chane it obra todo
                        {
                            found = true;
                            break;
                        }
                    }

                    if(found) nb_entry_correct_class ++;
                    else{

                        cout<<"\n found seq pbm: ("<<found_seq_start<<" , "<< found_seq_end<< ")"<<endl;
                        nb_entry_incorrect_class++;
                    }
                }
            }
            cout<<endl;

            std::cout << "inF : " << list_start_end_in_genome.at(i).first;
            std::cout << " , " << list_start_end_in_genome.at(i).second.size();
            std::cout << " || ";
            std::cout << "outF : " << list_res_parsedFile.at(i).first;
            std::cout << " , " << list_res_parsedFile.at(i).second.size();
            std::cout << "\n";
        }
        else // search it in the other bukket , i have to change this to map
        {

        }
    }


    cout<<" nb_entry_incorrect_class = "<<nb_entry_incorrect_class<<endl;
    cout<<" nb_entry_correct_class = "<<nb_entry_correct_class<<endl;

}

int main(int argc, char *argv[])
{
    //std::string file = R"(C:\Users\ibra\Desktop\Infernal\class_res_nbF_450_nbSeqs_min_3_max_4)";
    //std::string file = R"(C:\Users\ibra\Desktop\Infernal\class_res_nbF_450_nbSeqs_min_3_max_4_nohmm)";
    //std::string file = R"(C:\Users\ibra\Desktop\Infernal\class_res_nbF_450_nbSeqs_min_3_max_4_hmmonly)";

    //std::string file = R"(C:\Users\ibra\Desktop\Infernal\Clans ncRNA\class_res_Clans)";
    std::string file = R"(C:\Users\ibra\Desktop\Infernal\Clans ncRNA\class_res_Clans_nohmm)";
    //std::string file = R"(C:\Users\ibra\Desktop\Infernal\Clans ncRNA\class_res_Clans_hmmonly)";
    auto res = classParser(file);
    return 0;

    std::string file_to_parse = R"(C:\Users\ibra\Desktop\Infernal\nbF-10_nbSeqs-[30-35]\gS_noss_res_test_20)";
    std::string dir_test = R"(C:\Users\ibra\Desktop\Infernal\nbF-20_nbSeqs-[30-35]\Test)";

    //std::vector< std::pair<std::string, std::vector<std::pair<uint32_t , uint32_t >> > >

    check_search_result(dir_test, file_to_parse);
    return 0;
}