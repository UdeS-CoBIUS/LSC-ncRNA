//
// Created by ibra on 5/5/2020.
//

#ifndef MOTIFSEXTRACTIONSELECTION_COMMONMOTIFS_H
#define MOTIFSEXTRACTIONSELECTION_COMMONMOTIFS_H


#include <map>
#include <forward_list>
#include "SuffixTree_QuadraticTime.h"

using umap_seqId_nbOccs_ = unordered_map<unsigned int, unsigned int>;

class CommonMotifs {

    const SuffixTree_QuadraticTime &gst; // reference to Generalized Suffix Tree

    vector< umap_seqId_nbOccs_ > list_cm_umap_seqId_nbOccs; // contain all the umap for extracted motifs

    vector<string> list_cms; // list of all common Motifs

    vector<vector<uint32_t >> matrix_nbOcrrs_cmsSeqsIds; // M[i][j]=M[row][clm]=M[seqId][CMs]==nb_occurrences

    //std::map<unsigned int,unsigned int> map_str_len_idx; // map<str_lenght, str_idx>, the goal is to retrive the idx of all str, that are lower or upper bound to a given str lenght.
    vector<forward_list<uint32_t>> vec_flist; // vec_flist(max_len_str), each flist keep index of cm (of lenght i) in list_cm. the goal is to retrive the idx of all str, that are lower or upper bound to a given str lenght, to delete it if it is substring.

    vector<uint64_t > list_hashValue_cm_umap;

    //uint32_t nb_total_cm; // for debug
    //uint32_t nb_deleted_cm_not_accepted; // for debug
    //uint32_t nb_deleted_cm_from_list_cm; // for debug
    //uint32_t nb_unique_child; // for debug

    //----------------------- Selection parameters: --------------------------------------------------------------------
    unsigned int Beta=0; // BETA is percentage of finding cm in the family
    int Alpha=-1; // Alpha the variance in the number of occurences between seqs
    unsigned int nbOccrs_allowed=0;

    // -------  for debug and analysis
    unsigned int sum_cm_nb_occrs_alpha=0; // the value of repetition of an accpted motif against alpha


private: // methods

    umap_seqId_nbOccs_ cmsExtractionSelection_(const Node *node, const string& prefix);
    umap_seqId_nbOccs_ cmsExtractionSelectionDeletionCSMs_(const Node *node, const string& prefix);

    inline bool is_sub_motif(unsigned int idx_sm, unsigned int idx_m);
    inline bool is_sub_motif_hash(unsigned int idx_sm, unsigned int idx_m);
    inline bool is_sub_motif_hash(const string &ncm, const umap_seqId_nbOccs_ &ncm_umap, uint64_t ncm_hashValue, uint32_t idx_m);
    inline bool is_ncm_super_motif_hash(const string &ncm, const umap_seqId_nbOccs_ &ncm_umap, uint64_t ncm_hashValue, uint32_t idx_m);
    bool is_ncm_sub_motif_list_cms(const string &ncm, const umap_seqId_nbOccs_ &ncm_umap, uint64_t ncm_hashValue);

    void delete_sub_motifs_of_ncm_from_list_cms(string &ncm, umap_seqId_nbOccs_ &ncm_umap, unsigned long long int ncm_hashValue);
    inline void add_ncm_to_lists_cms(string &ncm, umap_seqId_nbOccs_ &ncm_umap, unsigned long long int ncm_hashValue);

    void update_valid_lists_cms_infos();

    void cleanUnnecessaryMemoryUsage(); // delete unnecessary intermediate memory usage.

    // ---------------------------------------------------------------------------------------------------------------------
    // ------------------------------------ CMs selection and filtrate Methods

    bool is_cm_accepted_according_To_selection_parameters(const string &cm, const umap_seqId_nbOccs_ &cm_umap_seqId_nbOccs);

    bool is_cm_accepted_only_for_one_family_beta_nbOcrrs(
            const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs,
            uint32_t percentage_same_family,
            uint32_t nbOccrs_allowed_local);

    bool
    is_cm_accepted_only_for_one_family_alpha_beta(
            const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs, uint32_t percentage_same_family,
            uint32_t alpha, uint32_t nbOccrs_allowed_local);

    // ---------------------------------------------------------------------------------------------------------------------
    // ------------------------------------ Utils Methods
    static void umapMerge(umap_seqId_nbOccs_ &des, umap_seqId_nbOccs_ &src);
    unordered_map<uint32_t , uint32_t > groupeCountSeqsByFamily(const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs);

public:
    CommonMotifs(const SuffixTree_QuadraticTime &gst, unsigned int beta);
    CommonMotifs(const SuffixTree_QuadraticTime &gst, unsigned int beta, int alpha);
    CommonMotifs(const SuffixTree_QuadraticTime &gst, unsigned int beta, int alpha, unsigned int nbOccrs_allowed);

    void cmsExtractionSelection();
    void cmsExtractionSelectionDeletionCSMs();

    void deleteSubMotifs_at_end_ES(); // delete sub-motifs after the phase of extraction and selection "cmsExtractionSelection()"
    void deleteSubMotifs_at_end_ES_usingHash(); // delete sub-motifs after the phase of extraction and selection "cmsExtractionSelection()" using Hashage

    void generateMatrixCmsSeqs();

    void print_list_cm_umap_seqId_nbOccs() const;

    const unsigned int get_nb_totals_cms() const;

    void print_matrix_cms_seqsId() const;

    void saveMatrixCMS_ToCsv_File(string file_output) const;
    void saveMatrixCMS_ToCsv_File_stdtostring(string file_output) const;
    void saveMatrixCMS_ToCsv_File_stdtostring_size_estimated(string file_output) const;
    void saveMatrixCMS_ToCsv_File_ofstream_write(string file_output) const;
    void saveMatrixCMS_ToCsv_File_dircrly(string file_output) const;
    void saveMatrixCMS_ToCsv_File_rowByrow(string file_output) const;
    void saveMatrixCMS_ToCsv_File_Chunked(string file_output, unsigned int chunk_size) const;

    ///bool can_delete_subMotif_hash(unsigned int idx_sm, unsigned int idx_lm);
    void
    flist_erase_all_idx_of_subMotif_of_ncm(string &ncm, umap_seqId_nbOccs_ &ncm_umap,
                                           unsigned long long int ncm_hashValue, forward_list<uint32_t> &flist);

    inline uint32_t get_nb_elements_in_vec_flist();

    void print_infos() const;
};


#endif //MOTIFSEXTRACTIONSELECTION_COMMONMOTIFS_H
