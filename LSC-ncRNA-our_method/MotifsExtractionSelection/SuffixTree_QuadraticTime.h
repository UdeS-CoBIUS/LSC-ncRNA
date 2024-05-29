//
// Created by ibra on 11/20/2019.
//

#ifndef SUFFIX_TREE_SUFFIXTREE_QUADRATICTIME_H
#define SUFFIX_TREE_SUFFIXTREE_QUADRATICTIME_H

#include <iostream>
#include <string>
#include <vector>
#include "Node.h"

using namespace std;

class SuffixTree_QuadraticTime {

private:
    Node *root_suffix_tree;

    // using either seaunecs without initial familes, or sequences are grouped in families
    vector<string> list_seqs; // 1) all sequences independent from initial clusturing
    vector<vector<string>> list_families_seqs;     // 2) sequnces are grouped in initial families (cluster)
    vector<unsigned int> list_sum_nb_seqs; // for indexing seqs by family

    bool sequences_are_grouped_in_families;

    bool use_index_for_seqs_id;

    unsigned int nb_all_sequences;

    uint32_t max_length_motif; // max length of motifs [1.. sequence_length], if 0, we use for each separate sequences its length.
    uint32_t min_length_motif;
public:
    uint32_t getMinLengthMotif() const;
    uint32_t getMaxLengthMotif() const;

private:
    // min length of motifs [1..(sequence_length-1)]. min_length < max_length

    bool use_min_max_length_motif;

private:
    void AddString(string str, unsigned int idx_seq);
    void addAllSubMotifMinMax(const string &str, uint32_t idx_str);
    void addAllSubMotifMin(const string &str, uint32_t idx_str);
    inline void addSubStr(const string &str, uint32_t idx_str, uint32_t from, uint32_t to);

    void set_Use_OrNot_IndexForSeqsId();
    void setIndexFamilies();

    pair<unsigned int, unsigned int> get_FamilyId_And_SeqId_SeqIdBased(unsigned int seq_id) const;
    pair<unsigned int, unsigned int> get_FamilyId_And_SeqId_IndexBased(const vector<unsigned int> &list_sum_nb_elements, unsigned int seq_id) const;

    unsigned int generateSeqId(unsigned int idx_family, unsigned int idx_seq_in_family) const;

    int get_num_char(char my_char) const;

    void help_initializer();


public:
    SuffixTree_QuadraticTime();

    SuffixTree_QuadraticTime(uint32_t maxLengthMotif, uint32_t minLengthMotif);


    virtual ~SuffixTree_QuadraticTime();

    const vector<string> &getListSeqs() const;
    const vector<vector<string>> &getListFamiliesSeqs() const;

    bool isSequencesAreGroupedByFamilies() const;
    bool isUseIndexForSeqsId() const;

    const string &getSeq(unsigned int id_seq) const;

    pair<unsigned int, unsigned int> get_FamilyId_And_SeqId(unsigned int seq_id) const;

    void GenerateSuffixTree(string &str);
    void GenerateGeneralizedSuffixTree(vector<string> &list_strings);
    void GenerateGeneralizedSuffixTree(vector<vector<string>> &list_grouped_strings);

    const Node *getRootSuffixTree() const;

    bool CheckSuffixTree(string str) const;
    bool CheckGST(vector<string> &listStrings) const;

    unsigned int getNbAllSequences() const;
};


#endif //SUFFIX_TREE_SUFFIXTREE_QUADRATICTIME_H
