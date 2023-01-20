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


void getMainArgv(int argc, char *argv[], string &dir_name, uint32_t &nb_families, uint32_t &min_nb_seqs_allowed,
                 uint32_t &max_nb_seqs_allowed, uint32_t &min_length_motif, uint32_t &max_length_motif,
                 uint32_t &is_delete_subMotifs, uint32_t &Beta)
{
    uint32_t nb_arg = 8;

    if(argc==nb_arg+1)
    {
        dir_name = argv[1];

        nb_families = strtol(argv[2], nullptr, 10);
        min_nb_seqs_allowed = strtol(argv[3], nullptr, 10);
        max_nb_seqs_allowed = strtol(argv[4], nullptr, 10);

        min_length_motif = strtol(argv[5], nullptr, 10);
        max_length_motif = strtol(argv[6], nullptr, 10);

        is_delete_subMotifs = strtol(argv[7], nullptr, 10);

        Beta = strtol(argv[8], nullptr, 10);
    }


    /*
    unordered_map<string,uint32_t > umap_main_argv;// text indice -d ---> position in argv
    umap_main_argv.reserve(6);
    umap_main_argv.insert(make_pair("-d",0));  // for dir_name
    umap_main_argv.insert(make_pair("-n",0));  // for nb files
    umap_main_argv.insert(make_pair("-mins",0)); // for min_nb_seqs_allowed
    umap_main_argv.insert(make_pair("-mans",0)); // for max_nb_seqs_allowed
    umap_main_argv.insert(make_pair("-milm",0)); // for min_length_motif
    umap_main_argv.insert(make_pair("-malm",0));  // for man_length_motif

    for (uint32_t i = 1; i < argc; ++i) { // i==0 is the programe name

        if(umap_main_argv.find(argv[i])==umap_main_argv.end())
        {
            cout<<"error in : arg ("<<(i+1)<<") {"<<argv[i]<<"} "<<endl;
        }
        else
        {
            umap_main_argv[argv]
        }
    }
    */


}

void print_test_vec(vector<forward_list<uint32_t>> &vec_flist, vector<string> &list_str)
{
    cout<<"\n --------------------- "<<endl;
    uint32_t len=0;
    for (auto &flist : vec_flist)
    {
        if(!flist.empty())
        {
            cout<<"len("<<len<<"): ";

            for (uint32_t i : flist)
            {
                cout<<"("<<list_str.at(i)<<", <"<<list_str.at(i).length()<<","<<i<<">) , ";
            }

            cout<<endl;
        }

        len++;
    }
    cout<<"\n --------------------- "<<endl;
}

void print_test_map(map<uint32_t , vector<uint32_t >> &map_str_len_listIdx, vector<string> &list_str)
{
    cout<<"\n --------------------- "<<endl;
    for (auto &it : map_str_len_listIdx)
    {
        cout<<"len("<<it.first<<"): ";

        for (uint32_t i : it.second)
        {
            cout<<"("<<list_str.at(i)<<", <"<<list_str.at(i).length()<<","<<i<<">) , ";
        }

        cout<<endl;
    }
    cout<<"\n --------------------- "<<endl;
}

bool remove_pred(uint32_t val)
{
    return val == 30;
}

void test_erase_from_forward_list()
{
    std::forward_list<int> mylist = {10, 10, 20, 30, 40, 10, 50, 10, 60};

    uint32_t skey = 10;

    for (int & val : mylist) std::cout<<val<<", ";

    auto prev = mylist.before_begin();
    for (auto it = mylist.begin(); it!=mylist.end(); )
    {
        if(*it==skey)
        {
            it = mylist.erase_after(prev);
            // break; // comment or uncomment to deal with only the first or all found element(s).
        }
        else
        {
            prev = it;
            ++it;
        }
    }

    std::cout<<"\n after deletion : "<<std::endl;
    for (int & val : mylist) cout<<val<<", ";
}

void flist_erase_all(forward_list<uint32_t> &flist,vector<string> &list_str, string str)
{

    auto prev = flist.before_begin();
    for (auto it = flist.begin(); it!=flist.end(); )
    {
        if(list_str.at(*it)==str)
        {
            it = flist.erase_after(prev);
            // break; // comment or uncomment to deal with only the first or all found element(s).
        }
        else
        {
            prev = it;
            ++it;
        }
    }
}

void test_map()
{

    vector<string> list_str{"ABCD","ABCDE","ibra","ibr","nasa","naceraty","fati","ibr","bily"};

    uint32_t min_ml = 2;
    uint32_t max_ml = 20;

    map<uint32_t , vector<uint32_t> > map_str_len_listIdx;
    vector<forward_list<uint32_t> > vec_strLen_FlistIdx(max_ml,forward_list<uint32_t>());

    uint32_t key;

    for (uint32_t i = 0; i < list_str.size(); ++i)
    {
        key= list_str.at(i).length();
        cout<<"("<<list_str.at(i)<<", <"<<list_str.at(i).length()<<","<<i<<">) , ";

        vec_strLen_FlistIdx.at(key).push_front(i);

        /*
        auto it = map_str_len_listIdx.find(key);
        if(it==map_str_len_listIdx.end()) // not found
        {
            map_str_len_listIdx.insert(make_pair(list_str.at(i).length(), vector<uint32_t >(1,i) ));
        }
        else // found
        {
            it->second.push_back(i);
        }
        */
    }

    cout<<"vec size = "<<list_str.size()<<endl;
    //cout<<"map size = "<<map_str_len_listIdx.size()<<endl;

    //print_test_map(map_str_len_listIdx,list_str);
    print_test_vec(vec_strLen_FlistIdx,list_str);

    // delete :
    // vector contain nb elements > 1  ==> delete that element
    // vector conatin nb elements = 1  ==> delete that element ++> vec empty
    //                                 ++> we can delete the entry from the map (pair(len,vec))
    //                                 +++> waste time in delete, re-arnge the map, and after that with higth propabilty we will insert the same key
    //                                 +++> don't delete the entry.
    //
    // vector conatin nb elements = 0 empty : do nothing, the key exist, but there is no index

    // as example we search str if it exist in map
    // 1) get list that have the same len
    // 2) compare with word by index from vec

    string str="ibr";
    flist_erase_all(vec_strLen_FlistIdx.at(str.length()),list_str,str);

    /*
    auto it = map_str_len_listIdx.find(str_len);
    if(it != map_str_len_listIdx.end())
    {
        for (unsigned int i=0; i < it->second.size(); ++i) // vector
        {
            unsigned int idx = it->second.at(i);

           if(str==list_str.at(idx))
           {
               it->second.erase(it->second.begin() + i);
               --i;
           }
        }
    }
    */

    cout<<"after deleting"<<endl;
    //print_test_map(map_str_len_listIdx,list_str);
    print_test_vec(vec_strLen_FlistIdx,list_str);

    //------------------------------
    flist_erase_all(vec_strLen_FlistIdx.at(str.length()),list_str,str);

    cout<<"after deleting"<<endl;
    //print_test_map(map_str_len_listIdx,list_str);
    print_test_vec(vec_strLen_FlistIdx,list_str);

    // insert at empty vector -------------------
    key= str.length();
    list_str.push_back(str);
    int idx_str = list_str.size()-1;

    vec_strLen_FlistIdx.at(key).push_front(idx_str);

    cout<<"after inserting"<<endl;
    //print_test_map(map_str_len_listIdx,list_str);
    print_test_vec(vec_strLen_FlistIdx,list_str);
}

void test_BST()
{
    std::vector<int> vec = { 10 , 115 ,15, 40, 35 , 55 , 75 };
    int searched_key = 40;

    std::set<int> my_BST(vec.begin(),vec.end());

    auto lower_bound = my_BST.lower_bound(searched_key);
    auto upper_bound = my_BST.upper_bound(searched_key);

    cout<<" Value Before = " << *(lower_bound.operator--())<<endl;
    cout<<" Value After = " << *upper_bound<<endl;
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