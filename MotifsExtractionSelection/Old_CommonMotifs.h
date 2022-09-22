//
// Created by ibra on 1/29/2020.
//

#ifndef ALLCOMMONMOTIFSNSTRINGS_COMMONMOTIFS_H
#define ALLCOMMONMOTIFSNSTRINGS_COMMONMOTIFS_H

#include <map>
#include "SuffixTree_QuadraticTime.h"

using umap_seqId_nbOccs_ = unordered_map<unsigned int, unsigned int>;

class Old_CommonMotifs {

    unsigned int Beta;

    const SuffixTree_QuadraticTime &gst;

    // we can add vector to store id familly, and or id seq in family.

    vector< umap_seqId_nbOccs_ > list_cm_umap_seqId_nbOccs; // just to help me and to debug, i should delete it after

    vector<vector<uint32_t >> matrix_commonMotifs_seqsId; // we don't need the size, because the vector give us the size.
    //seqId[Old_CommonMotifs[nb_occurrences]]

    // list motifs: list_of_all_commonMotifs
    //  motifs ids-->
    // |   [][j][j][j][j][j][]
    // |   [i]{k}
    // V   [i]
    //     [i]
    // seqId
    //
    // i: les sequences
    // j: commons motifs
    // k:[i,j] : a vector of position of the common motif(j) in the sequence i.



    std::map<unsigned int,unsigned int> map_str_len_idx; // map of <str_lenght, str_idx>,
    // the goal is to retrive the idx of all str,
    // that are lower or upper bound to a given str lenght.

    std::vector<bool > bitVec_motifs; // keep track if a motif is deleted or not

    umap_seqId_nbOccs_ fillMatrixCommonMotifsSeqsId(const Node *node, const string& prefixe);
    umap_seqId_nbOccs_ fillMatrixCommonMotifsSeqsId_deleteSubMotifs(const Node *node, const string& prefixe);

    static void umapMerge(unordered_map<unsigned int, unsigned int> &des,
                          unordered_map<unsigned int, unsigned int> &src);

    void insert_cm_umap_seqId_nbOccs_In_matrix_commonMotifs_seqsId(
            const unordered_map<unsigned int, unsigned int> &umap_seqId_nbOccs);

    void insert_cm_umap_seqId_nbOccs_In_matrix_commonMotifs_seqsId_deleteSubMotifs(
            const unordered_map<unsigned int, unsigned int> &umap_seqId_nbOccs);

    bool isListSeqIdAreInSameFamily(const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs);

    bool
    isCm_Common_For_ListSeqId_Percentage(const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs,
                                         uint32_t percentage_same_family,
                                         uint32_t percentage_multiple_families);

    bool isCm_Common_For_ListSeqId_Percentage(const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs,
                                              uint32_t percentage_same_family,
                                              uint32_t percentage_multiple_families,
                                              uint32_t alpha);

public:

    Old_CommonMotifs(const SuffixTree_QuadraticTime &gst, unsigned int beta);

    void cleanUnnecessaryMemoryUsage(); // delete unnecessary intermediate memory usage.

    void generateMatrixCmsSeqs();
    void generateMatrixCmsSeqs_deleteSubMotifs();

   void print_list_cm_umap_seqId_nbOccs() const;

    void print_matrix_commonMotifs_seqsId() const;

    void saveToCsvMatrixCommonMotifsSeqsIds();

    void saveToCsvMatrixCommonMotifsSeqsIds_buffer_reserve();

    void saveToCsvMatrixCommonMotifsSeqsIds_buffer();

    bool
    isCm_Common_For_ListSeqId_Percentage_cmInter(
            const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs, uint32_t percentage_same_family);

    unordered_map<uint32_t, uint32_t>
    groupeCountSeqsByFamily(const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs);

    bool
    isCm_Common_For_ListSeqId_Percentage_cmInter_2(
            const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs,
            uint32_t percentage_same_family);

    bool areCommonMotifsInSameSeqs(
            const unordered_map<unsigned int, unsigned int> &left_cm_umap,
            const unordered_map<unsigned int, unsigned int> &right_cm_umap);

    bool areCommonMotifsInSameSeqs(uint32_t lcm, uint32_t rcm);

    void deletSubMotifsFromTheMatrix();

    bool compare(uint32_t i, uint32_t j);

    bool can_delete_subMotif(unsigned int idx_sm, unsigned int idx_m);



    void delete_upDate_SubMotifsFromTheMatrix();

    vector<string> list_of_all_commonMotifs;

    bool is_cm_sub_motifs(const string& new_cm, const umap_seqId_nbOccs_& umap_new_cm);

    bool areCommonMotifsInSameSeqs_sameNbOccrs(uint32_t lcm, uint32_t rcm);
    bool areCommonMotifsInSameSeqs_sameNbOccrs(umap_seqId_nbOccs_ umap_new_cm, uint32_t idx_cm_mtrx);

    void delete_subMotifs(const string &new_cm, const umap_seqId_nbOccs_ &umap_new_cm);
};


#endif //ALLCOMMONMOTIFSNSTRINGS_COMMONMOTIFS_H
