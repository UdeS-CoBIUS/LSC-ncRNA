//
// Created by ibra on 1/29/2020.
//

#include <fstream>
#include <unordered_set>
#include <numeric>
#include <algorithm>
#include "Old_CommonMotifs.h"
#include "SubSetDistancePercentage.h"

void Old_CommonMotifs::generateMatrixCmsSeqs()
{
    string root_motif; // at the root, we have an empty motif.

    //fillMatrixCommonMotifsSeqsId(this->gst.getRootSuffixTree(), root_motif);

    const Node * root = this->gst.getRootSuffixTree();

    for (auto & child : root->listChildren)
    {
        if(child)
        {
            fillMatrixCommonMotifsSeqsId(child, root_motif);
        }
    }
}

umap_seqId_nbOccs_ Old_CommonMotifs::fillMatrixCommonMotifsSeqsId(const Node *node, const string& prefixe)
{
    string motif_active_node = prefixe + node->cara;
    umap_seqId_nbOccs_ umap_seqId_nbOccs_active_node{};

    //if(motif_active_node.size() >= this->gst.getMinLengthMotif())
    if(motif_active_node.size() > this->gst.getMinLengthMotif())
    {
        if(!node->map_seqId_nbOccs.empty())
            umap_seqId_nbOccs_active_node = node->map_seqId_nbOccs; // COPY

        umap_seqId_nbOccs_ umap_seqId_nbOccs_sub_tree;

        for (auto & child : node->listChildren)
        {
            if(child)
            {
                umap_seqId_nbOccs_sub_tree = fillMatrixCommonMotifsSeqsId(child, motif_active_node);
                umapMerge(umap_seqId_nbOccs_active_node, umap_seqId_nbOccs_sub_tree);
            }
        }

        // at this stage, we have the umap_seqId_nbOccs_active_node that contain all seqId and position of the actuel node,
        // i.e: where the motifs that end in this node
        // add this common motifs to the matrix.

        //if(umap_seqId_nbOccs_active_node.size()>GAMMA)
        //if(!umap_seqId_nbOccs_active_node.empty() && isCm_Common_For_ListSeqId_Percentage(umap_seqId_nbOccs_active_node,BETA,BETA_2))
        //if(!umap_seqId_nbOccs_active_node.empty() && isCm_Common_For_ListSeqId_Percentage_cmInter(umap_seqId_nbOccs_active_node,BETA))
        if(!umap_seqId_nbOccs_active_node.empty() && isCm_Common_For_ListSeqId_Percentage_cmInter_2(umap_seqId_nbOccs_active_node,Beta))
        {
            list_of_all_commonMotifs.push_back(motif_active_node);

            insert_cm_umap_seqId_nbOccs_In_matrix_commonMotifs_seqsId(umap_seqId_nbOccs_active_node);

            // important:
            // this code is important , not for here, but whene  i will construct the Htree directly without distMatrix, so i will need only the seqs that share the CM.
            //list_cm_umap_seqId_nbOccs.push_back(umap_seqId_nbOccs_active_node);
        }
    }
    else
    {
        for (auto & child : node->listChildren)
        {
            if(child)
            {
                fillMatrixCommonMotifsSeqsId(child, motif_active_node);
            }
        }
    }

    return umap_seqId_nbOccs_active_node;
}

void Old_CommonMotifs::umapMerge(umap_seqId_nbOccs_ &des, umap_seqId_nbOccs_ &src)
{
    for (auto &element:src)
    {
        if(des.find(element.first)==des.end())
        {
            des.insert(element);
        }
        else
        {
            // let the word 0:ABCABC, (we consider ABC as one):
            //          root
            //           |
            //          ABC
            //         /  |
            //    {0,1})  |
            //           ABC({0,1})
            // that mean ABCABC exist 1 times
            // the motifs ABC exist 2 times, because it is itself leaf and father (prefix) of ABC.
            // as result:
            // root-->ABC({0,2})-->ABC({0,1})

            des[element.first]+=element.second; // add the nb occs of its children to him
        }
    }
}

void Old_CommonMotifs::insert_cm_umap_seqId_nbOccs_In_matrix_commonMotifs_seqsId(
        const umap_seqId_nbOccs_ &umap_seqId_nbOccs)
{
    vector<uint32_t > list_seqId_nbOccs(this->gst.getNbAllSequences());

    for(const auto &pair_seqId_nbOccs:umap_seqId_nbOccs)
    {
        uint32_t i = pair_seqId_nbOccs.first;
        uint32_t val = pair_seqId_nbOccs.second;

        list_seqId_nbOccs.at(i) = val;
    }

    this->matrix_commonMotifs_seqsId.push_back(list_seqId_nbOccs);
}



void Old_CommonMotifs::print_list_cm_umap_seqId_nbOccs() const
{
    for (uint32_t i=0 ; i<list_cm_umap_seqId_nbOccs.size(); ++i)
    {
        const auto & motif = list_of_all_commonMotifs.at(i);
        const auto & cm_umap_seqId_nbOccs = list_cm_umap_seqId_nbOccs.at(i);

        //1) get common motif of the map:
        cout<<motif;
        cout<<" ==> ";

        //2) print for each paire: seqId and nbOccs:
        for (const auto & pair_seqId_nbOccs:cm_umap_seqId_nbOccs)
        {
            cout<<"[ ";
            cout << pair_seqId_nbOccs.first; // seqId
            cout<<" : ";
            cout<<util::to_string(pair_seqId_nbOccs.second); // nbOccs
            cout<<" ]";
        }

        cout<<endl;
    }
}

void Old_CommonMotifs::print_matrix_commonMotifs_seqsId() const
{
    unsigned int nb_seqs= this->gst.getNbAllSequences();
    unsigned int nb_motifs = this->list_of_all_commonMotifs.size();

    // 1) print the motifs
    cout<<endl<<"        ";
    for (const auto & motif:this->list_of_all_commonMotifs)
    {
        cout<<motif <<"|  ";
    }
    cout<<endl;

    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        cout<<"seq("<<i<<"): ";

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            cout<<this->matrix_commonMotifs_seqsId[j][i];
            cout<<"| ";
        }

        cout<<"\n____________"<<endl;
    }


}

void Old_CommonMotifs::saveToCsvMatrixCommonMotifsSeqsIds()
{
    unsigned int nb_seqs= this->gst.getNbAllSequences();
    unsigned int nb_motifs = this->list_of_all_commonMotifs.size();

    ofstream file_csv;
    string csv_separator =",";
    file_csv.open ("matrix_motifs_seqIds.csv");
    cout<<" file_csv.open (\"matrix_motifs_seqIds.csv\"); "<<endl;

    // 1) index is the global seqId, seqIdInFam is index of a given seq in its family, and familyId is the class.
    file_csv<<"index"<<csv_separator
            <<"seqIdInFam"<<csv_separator
            <<"familyId"<<csv_separator;

    // 2) print all the motifs in one row.
    for (const auto & motif:this->list_of_all_commonMotifs)
    {
        file_csv<<motif
                <<csv_separator;
    }
    file_csv<<endl;


    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->gst.isSequencesAreGroupedByFamilies())
        {
            auto pair_familyId_SeqId = this->gst.get_FamilyId_And_SeqId(i); //(idx_family,idx_seq_in_family)

            // global seq id , the id are soposed to be index based , meaning all sequnces in all families have the same continues system of id, i,i+1...
            file_csv<<i<<csv_separator
                    <<pair_familyId_SeqId.second<<csv_separator //idx_seq_in_family
                    <<pair_familyId_SeqId.first<<csv_separator; //idx_family
        } else
        {
            file_csv<<i<<csv_separator
                    <<i<<csv_separator //idx_seq_in_family is the same as global index
                    <<0<<csv_separator; //idx_family is 0, we have one family
        }

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            file_csv<<matrix_commonMotifs_seqsId[j][i]; // nb occurence
            file_csv<<csv_separator;
        }
        file_csv<<endl;
    }

    file_csv.close();

    cout<<"file_csv.close();"<<endl;
}


void Old_CommonMotifs::saveToCsvMatrixCommonMotifsSeqsIds_buffer_reserve()
{
    unsigned int nb_seqs= this->gst.getNbAllSequences();
    unsigned int nb_motifs = this->list_of_all_commonMotifs.size();

    cout<<" Total nb_seqs = "<<nb_seqs<<endl;
    cout<<" Total nb_motifs = "<<nb_motifs<<endl;

    ofstream file_csv;
    const string csv_separator =",";
    string buffer;
    uint32_t size_buffer = 0;

    size_buffer += nb_motifs * nb_seqs; // the size of the matrix_commonMotifs_seqsId

    // this supose that we have only one digit, but we can have tow digit!! so howe  i should do this
    // generaly we have only one digit, but tow digit exist.
    // i think to add 1/3 of the space to cover the possible tow digit
    size_buffer += size_buffer/3;


    size_buffer += 3*nb_seqs; // we add 3 colomn {index, seqIdInFam, familyId}
    size_buffer += nb_seqs; // each seq is a row, so it is an "\n" at the end

    // 1) index is the global seqId, seqIdInFam is index of a given seq in its family, and familyId is the class.
    string idx_colo_name;
        idx_colo_name+="index";
        idx_colo_name+=csv_separator;
        idx_colo_name+="seqIdInFam";
        idx_colo_name+=csv_separator;
        idx_colo_name+="familyId";
        idx_colo_name+=csv_separator;


    uint32_t size_all_colomn_name = idx_colo_name.size();


    // 2) compute the total size of all motifs
    for (const auto & motif:this->list_of_all_commonMotifs)
    {
        size_all_colomn_name+=motif.size()+1; // +1 for the separtor
    }
    size_all_colomn_name +=1; // for the endline.

    size_buffer+=size_all_colomn_name;

    buffer.reserve(size_buffer);

    buffer+=idx_colo_name;
    // 2) compute the total size of all motifs
    for (const auto & motif:this->list_of_all_commonMotifs)
    {
        buffer+=motif;
        buffer+=csv_separator;
    }
    buffer+="\n";

    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->gst.isSequencesAreGroupedByFamilies())
        {
            auto pair_familyId_SeqId = this->gst.get_FamilyId_And_SeqId(i); //(idx_family,idx_seq_in_family)

            // global seq id , the id are soposed to be index based , meaning all sequnces in all families have the same continues system of id, i,i+1...
            buffer+=util::to_string(i);
            buffer+=csv_separator;
            buffer+=util::to_string(pair_familyId_SeqId.second);
            buffer+=csv_separator; //idx_seq_in_family
            buffer+=util::to_string(pair_familyId_SeqId.first);
            buffer+=csv_separator; //idx_family
        } else
        {
            buffer+=util::to_string(i);
            buffer+=csv_separator;
            buffer+=util::to_string(i);
            buffer+=csv_separator; //idx_seq_in_family is the same as global index
            buffer+="0";
            buffer+=csv_separator; //idx_family is 0, we have one family
        }

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            buffer+=util::to_string(matrix_commonMotifs_seqsId[j][i]); // nb occurence
            buffer+=csv_separator;
        }
        buffer+="\n";
    }

    //-------------------------------------------------------------------
    //string dir_path = R"(C:\Users\ibra\OneDrive - USherbrooke\DatatSet\pr-2_nbf-X_MinNbSeq-2_MaxNbSeq-60_MinLm-2_MaxLm-10\Mtd-2_60prec_nbf-X\matrix_motifs_seqIds.csv)";
    string dir_path = R"(matrix_motifs_seqIds.csv)";
    //file_csv.open ("matrix_motifs_seqIds.csv");
    file_csv.open (dir_path);
    cout<<" file_csv.open (\"matrix_motifs_seqIds.csv\"); "<<endl;

    file_csv<<buffer;

    file_csv.close();
    cout<<"file_csv.close();"<<endl;
}

void Old_CommonMotifs::saveToCsvMatrixCommonMotifsSeqsIds_buffer()
{
    unsigned int nb_seqs= this->gst.getNbAllSequences();
    unsigned int nb_motifs = this->list_of_all_commonMotifs.size();

    ofstream file_csv;
    const string csv_separator =",";
    string buffer;

    // 1) index is the global seqId, seqIdInFam is index of a given seq in its family, and familyId is the class.
    buffer+="index";
    buffer+=csv_separator;
    buffer+="seqIdInFam";
    buffer+=csv_separator;
    buffer+="familyId";
    buffer+=csv_separator;


    // 2) compute the total size of all motifs
    for (const auto & motif:this->list_of_all_commonMotifs)
    {
        buffer+=motif;
        buffer+=csv_separator;
    }
    buffer+="\n";

    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->gst.isSequencesAreGroupedByFamilies())
        {
            auto pair_familyId_SeqId = this->gst.get_FamilyId_And_SeqId(i); //(idx_family,idx_seq_in_family)

            // global seq id , the id are soposed to be index based , meaning all sequnces in all families have the same continues system of id, i,i+1...
            buffer+=i;
            buffer+=csv_separator;
            buffer+=pair_familyId_SeqId.second;
            buffer+=csv_separator; //idx_seq_in_family
            buffer+=pair_familyId_SeqId.first;
            buffer+=csv_separator; //idx_family
        } else
        {
            buffer+=i;
            buffer+=csv_separator;
            buffer+=i;
            buffer+=csv_separator; //idx_seq_in_family is the same as global index
            buffer+="0";
            buffer+=csv_separator; //idx_family is 0, we have one family
        }

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            buffer+=matrix_commonMotifs_seqsId[j][i]; // nb occurence
            buffer+=csv_separator;
        }
        buffer+="\n";
    }

    //-------------------------------------------------------------------
    file_csv.open ("matrix_motifs_seqIds.csv");
    cout<<" file_csv.open (\"matrix_motifs_seqIds.csv\"); "<<endl;

    file_csv<<buffer;

    file_csv.close();
    cout<<"file_csv.close();"<<endl;
}

Old_CommonMotifs::Old_CommonMotifs(const SuffixTree_QuadraticTime &gst, unsigned int beta) : gst(gst) { this->Beta=beta;}

bool Old_CommonMotifs::isCm_Common_For_ListSeqId_Percentage(
        const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs,
        uint32_t percentage_same_family, uint32_t percentage_multiple_families)
{
    uint32_t nb_all_seqs_comparedTo = 0;
    uint32_t nb_seqs_have_cm = cm_umap_seqId_nbOccs.size();  // nb seqs that have the common motif.

    if(this->gst.isSequencesAreGroupedByFamilies())
    {
        if(isListSeqIdAreInSameFamily(cm_umap_seqId_nbOccs))
        {
            uint32_t global_seqId = cm_umap_seqId_nbOccs.begin()->first;
            uint32_t first_family_id = this->gst.get_FamilyId_And_SeqId(global_seqId).first; //(idx_family,idx_seq_in_family)

            nb_all_seqs_comparedTo = this->gst.getListFamiliesSeqs().at(first_family_id).size();

            return ( (double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo >= percentage_same_family );
        } else
        {
            // 1) get all the family id
            unordered_set<uint32_t > set_families_ids;
            uint32_t global_seqId;
            uint32_t family_id;

            for (const auto &cm_pair_seqId_nbOcc : cm_umap_seqId_nbOccs)
            {
                global_seqId = cm_pair_seqId_nbOcc.first;
                family_id = this->gst.get_FamilyId_And_SeqId(global_seqId).first; //(idx_family,idx_seq_in_family)

                if(set_families_ids.find(family_id)==set_families_ids.end())
                {
                    set_families_ids.insert(family_id);
                }
            }

            // 2) computes all size , and moyenne
            uint32_t moyenne=0;

            for (uint32_t familyId:set_families_ids)
            {
                nb_all_seqs_comparedTo += this->gst.getListFamiliesSeqs().at(familyId).size();
            }

            moyenne = nb_all_seqs_comparedTo / set_families_ids.size();

            // based all seqs in all families that share the cm
            //return ( (double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo >= percentage_multiple_families );

            // based on moyenne of all seqs in all families that share the cm
            return ( (double)(nb_seqs_have_cm*100)/(double)moyenne >= percentage_multiple_families );
        }
    }

    // here, meaning sequnces are not grouped, so they are supposed tu be all in the same family
    nb_all_seqs_comparedTo = this->gst.getNbAllSequences();
    return ( (double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo >= percentage_same_family );
}


/**
 *  a cm motif is reprenstatif of groupe of families if it is represntatif of each family alone according to BETA_2
 *  BETA_2 my be just BETA wich make sense (BETA is the thershould that say if cm is reprensttaif of one family)
 *    it make sense that BETA_2==BETA, becauese, simly if we veiw the thing like:
 *    let A, and B tow sets,
 *    let cm(A)={a1,a2,...,an} cm motif for the set A with each a_i is represntatif of A (mean exist on >=BETA of sequnces of A)
 *    let cm(B)={b1,b2,...,bm} cm motif for the set B with each b_i is represntatif of B (mean exist on >=BETA of sequnces of B)
 *    C=cm(A) \inter c(B) = {c1,c2,...,cp}
 *    this mean each cm ci is representativ of A and B in the same time
 *    when we construct the HTree, so the node father of A and B will containe only cm , i.e: C
 * @param cm_umap_seqId_nbOccs
 * @param percentage_same_family
 * @return
 */
bool Old_CommonMotifs::isCm_Common_For_ListSeqId_Percentage_cmInter(
        const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs, uint32_t percentage_same_family)
{
    uint32_t nb_all_seqs_comparedTo = 0;
    uint32_t nb_seqs_have_cm = cm_umap_seqId_nbOccs.size();  // nb seqs that have the common motif.

    if(this->gst.isSequencesAreGroupedByFamilies())
    {
        const auto umap_family_nbSeqs = groupeCountSeqsByFamily(cm_umap_seqId_nbOccs);

        for(const auto & pair_familyId_nbSeqs:umap_family_nbSeqs)
        {
            uint32_t first_family_id = pair_familyId_nbSeqs.first;
            nb_all_seqs_comparedTo = this->gst.getListFamiliesSeqs().at(first_family_id).size();

            if( ((double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo) < percentage_same_family )
                return false;
        }

        return true;
    }

    // here, meaning sequnces are not grouped, so they are supposed tu be all in the same family
    nb_all_seqs_comparedTo = this->gst.getNbAllSequences();
    return ( (double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo >= percentage_same_family );
}

/**
 *  for multiple family, if cm is stisfy only one family, we take it.
 *  because: the motifs are common between many sequnces at the same time,
 *  for example cm=ACTG = {1,2,3,4,5, 100,101} ,
 *  we supose that cm is represntatif of the family where the seqs {1,2,3,4,5} belong
 *  but is not represntative for the family where the seqs {100,101} belong.
 *  for that, if we refuse the motif,
 *  that mean we will miss the motif that represnt the first family
 *  and the pbm is the is no way to get it
 *  (unless we treat each family alone, to get cm motifs for it, and then to construct HCTee, we do intersection beween list of cm)
 *  so for that we accepte it, because the goal is to put it in the generated matrix for classification by Nabil (learning methodes)
 * @param cm_umap_seqId_nbOccs
 * @param percentage_same_family
 * @return
 */
bool Old_CommonMotifs::isCm_Common_For_ListSeqId_Percentage_cmInter_2(
        const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs,
        uint32_t percentage_same_family)
{
    uint32_t nb_all_seqs_comparedTo = 0;
    uint32_t nb_seqs_have_cm = cm_umap_seqId_nbOccs.size();  // nb seqs that have the common motif.

    if(this->gst.isSequencesAreGroupedByFamilies())
    {
        const auto umap_family_nbSeqs = groupeCountSeqsByFamily(cm_umap_seqId_nbOccs);

        for(const auto & pair_familyId_nbSeqs:umap_family_nbSeqs)
        {
            uint32_t first_family_id = pair_familyId_nbSeqs.first;
            nb_all_seqs_comparedTo = this->gst.getListFamiliesSeqs().at(first_family_id).size();

            if( ((double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo) >= percentage_same_family )
                return true;
        }

        return false;
    }

    // here, meaning sequnces are not grouped, so they are supposed tu be all in the same family
    nb_all_seqs_comparedTo = this->gst.getNbAllSequences();
    return ( (double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo >= percentage_same_family );
}


/**
 *  whene calling this method, we supose (and check before call) that sequnces are grouped in fiffirent families
 * @param cm_umap_seqId_nbOccs
 * @return
 */
bool Old_CommonMotifs::isListSeqIdAreInSameFamily(const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs)
{
    if(cm_umap_seqId_nbOccs.empty())
        return false;

    uint32_t global_seqId = cm_umap_seqId_nbOccs.begin()->first;
    uint32_t first_family_id = this->gst.get_FamilyId_And_SeqId(global_seqId).first; //(idx_family,idx_seq_in_family)

    uint32_t family_id;
    for (const auto &cm_pair_seqId_nbOcc : cm_umap_seqId_nbOccs)
    {
        global_seqId = cm_pair_seqId_nbOcc.first;
        family_id = this->gst.get_FamilyId_And_SeqId(global_seqId).first; //(idx_family,idx_seq_in_family)

        if(family_id != first_family_id) return false;
    }


    return true;
}


/**
 * check if common motifs left and common motif right are in the same sequences.
 * @param left_cm_umap left cm_umap_seqId_nbOccs : unordered_map of (seq id, nb occrs of common motif)  for a given cm
 * @param right_cm_umap rigth cm_umap_seqId_nbOccs : unordered_map of (seq id, nb occrs of common motif)  for a given cm
 * @return if the two list contain the same keys ==> if contain the same seqs ids.
 */
bool Old_CommonMotifs::areCommonMotifsInSameSeqs(
        const unordered_map<unsigned int, unsigned int> &left_cm_umap,
        const unordered_map<unsigned int, unsigned int> &right_cm_umap)
{

    return
    (left_cm_umap.size() != right_cm_umap.size())
    &&
    std::equal(left_cm_umap.begin(), left_cm_umap.end(), right_cm_umap.begin(),
                  [] (auto a, auto b) { return a.first == b.first; });
}

/**
 * check if lcm and rcm if they share the same seqs
 * the check is done in the matrix that hold the number of occurences
 * lcm and rcm share the same sequnces if thier seqs hold a nboccrs >0,
 * because, 0 mean that the cm is not exist in the seq.
 * or if they share the same seqs that hold 0.
 * @param lcm idx of the left cm in matrix (seqs id, motif) that hold nbOccrs
 * @param rcm idx of the right cm in matrix (seqs id, motif) that hold nbOccrs
 * @return if the lcm and rcm are in the same seqs
 */
bool Old_CommonMotifs::areCommonMotifsInSameSeqs(uint32_t lcm, uint32_t rcm)
{
    unsigned int nb_seqs= this->gst.getNbAllSequences();
    unsigned int nb_motifs = this->list_of_all_commonMotifs.size();

    // we don't have this case, else it is a pbm
    // 1) idx of lcm and rcm must <nb_motifs
    //if(lcm>=nb_motifs || rcm>=nb_motifs)  // pbm of index
    //    return false;

    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->matrix_commonMotifs_seqsId[lcm][i]==0
        &&
        this->matrix_commonMotifs_seqsId[rcm][i]!=0)
            return false;

        if(this->matrix_commonMotifs_seqsId[lcm][i]!=0
           &&
           this->matrix_commonMotifs_seqsId[rcm][i]==0)
            return false;
    }

    return true;
}

/**
 * check if lcm and rcm if they share the same seqs and the same nb of occrs
 * the check is done in the matrix that hold the number of occurences
 * lcm and rcm share the same sequnces and same nboccrs if thier seqs hold the same nboccrs
 * this condition is sufficant to ensure that are they in the same seqs,
 * because if they have no nboccrs, they must 0 wich mean they don't exist at the same time.
 * @param lcm idx of the left cm in matrix (seqs id, motif) that hold nbOccrs
 * @param rcm idx of the right cm in matrix (seqs id, motif) that hold nbOccrs
 * @return if the lcm and rcm are in the same seqs
 */
bool Old_CommonMotifs::areCommonMotifsInSameSeqs_sameNbOccrs(uint32_t lcm, uint32_t rcm)
{
    unsigned int nb_seqs= this->gst.getNbAllSequences();
    unsigned int nb_motifs = this->list_of_all_commonMotifs.size();

    // we don't have this case, else it is a pbm
    // 1) idx of lcm and rcm must <nb_motifs
    //if(lcm>=nb_motifs || rcm>=nb_motifs)  // pbm of index
    //    return false;

    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->matrix_commonMotifs_seqsId[lcm][i]
            != this->matrix_commonMotifs_seqsId[rcm][i])
            return false;
    }

    return true;
}

/**
 * to check if new_cm have the same sequnces and nb occrs
 * 1) we must find the seqsId (and numbers occrs) of map_new_cm in the clomn of the motifs
 * 2) the other case of the motifs must be 0
 * @param umap_new_cm  map of the new found cm
 * @param idx_cm_mtrx idx of coloumn in the matrix of the motifs that we compar with
 * @return if they have same sequnces and same number of occurences.
 */
bool Old_CommonMotifs::areCommonMotifsInSameSeqs_sameNbOccrs(umap_seqId_nbOccs_ umap_new_cm, uint32_t idx_cm_mtrx)
{
    unsigned int nb_seqs= this->gst.getNbAllSequences();

    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->matrix_commonMotifs_seqsId[idx_cm_mtrx][i]
           != 0)
        {
            if(umap_new_cm.find(i)!=umap_new_cm.end()) // the seqId exist
            {
                if( this->matrix_commonMotifs_seqsId[idx_cm_mtrx][i]
                != umap_new_cm[i])
                    return false;

            } else // the seq have nbOccrs (!=0) && it is not exist in new_cm
                return false;
        }
    }

    return true;
}

bool Old_CommonMotifs::can_delete_subMotif(unsigned int idx_sm, unsigned int idx_m)
{
    //1 check if the sm is submotif of m
    if(list_of_all_commonMotifs.at(idx_m).find(list_of_all_commonMotifs.at(idx_sm)) ==std::string::npos )
        return false;

    // at this stage, the sm is submotif of m.
    //2) check if they share the same sequnces.
    //return areCommonMotifsInSameSeqs(idx_sm,idx_m);

    // 3) check the number of occurences
    //return true;

    return areCommonMotifsInSameSeqs_sameNbOccrs(idx_sm,idx_m);
}

/**
 * delete the submotif from the matrix and the list of motifs as well
 * we delete if: 1) the sm and m share the same sequences and
 *               2) have the same number of occurences
 */
void Old_CommonMotifs::deletSubMotifsFromTheMatrix()
{
    unsigned int nb_deleted_elements=0;

    // 1) creat an index of all motifs, and sort it according to the motis lenght
    vector<int > list_motif_id(list_of_all_commonMotifs.size());
    std::iota(list_motif_id.begin(),list_motif_id.end(),0);

    std::sort(list_motif_id.begin(),list_motif_id.end(),[&](uint32_t i, uint32_t j){
        return list_of_all_commonMotifs.at(i).size() < list_of_all_commonMotifs.at(j).size();
    });

    // only check for debug
    //cout<<"before deleting : "<<endl;
    //for(unsigned int i =0; i<list_motif_id.size(); ++i) cout<<"| "<<list_of_all_commonMotifs.at(list_motif_id.at(i));

    // 2) check if motis i is sub-motifs of one of motifs that are greater than it. if yes delete it.
    for(unsigned int i =0; i<list_motif_id.size()-1;++i)
    {
        for (unsigned int j = i+1; j < list_motif_id.size(); ++j)
        {
            //if(areCommonMotifsInSameSeqs(list_motif_id.at(i),list_motif_id.at(j)))
            if(can_delete_subMotif(list_motif_id.at(i),list_motif_id.at(j)))
            //if(areCommonMotifsInSameSeqs_sameNbOccrs(list_motif_id.at(i),list_motif_id.at(j)))
            {
                // delete the element form list_common_motifs and the matrix at the idx i.
               //cout<<"["<<i<<", "<<j<<"] same"<<endl;

               // this dos't work, becease with deletion we have idx pbm...
               //list_of_all_commonMotifs.erase(list_of_all_commonMotifs.begin()+list_motif_id.at(i));

               list_motif_id.at(i)=-1; // to say that we don't need it.
               nb_deleted_elements++;

                break;

            }
        }
    }

//    cout<<"\nAfter deleting : "<<endl;
//    for(int i : list_motif_id){
//        if(i!=-1){
//            cout<<"| "<<list_of_all_commonMotifs.at(i);
//        }
//    }

    //3 ) update the matrix : keep only motifs that are not deleted in the matrix
    vector<string> tmp_list_of_all_commonMotifs;
    tmp_list_of_all_commonMotifs.reserve(list_of_all_commonMotifs.size()-nb_deleted_elements);

    vector<vector<uint32_t >> tmp_matrix_commonMotifs_seqsId;
    tmp_matrix_commonMotifs_seqsId.reserve(matrix_commonMotifs_seqsId.size()-nb_deleted_elements);


    for(int i : list_motif_id){
        if(i!=-1){
            tmp_list_of_all_commonMotifs.push_back( list_of_all_commonMotifs.at(i) );
            tmp_matrix_commonMotifs_seqsId.push_back( matrix_commonMotifs_seqsId.at(i) );
        }
    }

    //cout<<"\n tmp-list_strings : "<<endl;
    //cout<<util::to_string(tmp_list_of_all_commonMotifs);

    this->list_of_all_commonMotifs = tmp_list_of_all_commonMotifs;
    this->matrix_commonMotifs_seqsId = tmp_matrix_commonMotifs_seqsId;

    //this->print_matrix_commonMotifs_seqsId();

    cout<<"nb motifs after deletion : "<<list_of_all_commonMotifs.size()<<endl;
}

/**
 * delete the submotif from the matrix and the list of motifs as well
 * we delete if: 1) the sm and m share the same sequences and
 *               2) have the same number of occurences
 *
 * up date the number of occurences of a submotifs if:
 *  1) the sm and m share the same sequences and
 *  2) they don't have the same number of occurences
 *
 *  the firs idea was:
 *  -check if sm and m are in the same sequneces
 *  - if they have the same nb occrs --> delete
 *  - else
 *  - update nb occrs de sm by : nbOccrs(sm) = nbOccrs(sm) - nbOccrs(m)
 *   which seem to be logic, and this is true if there is only one motif that include that subm
 *   and it's remain true whith multiple motifs that are separate in diffrent position
 *   in this case we can do the same thing, substruct the diffirence.
 *   example, let sm, m1, m2 where sm is subm of m1 and also of m2, m1 and m2 are not overlaping or one include the other.
 *   in this case, nbOccrs(sm) = nbOccrs(sm) - nbOccrs(m1) - nbOccrs(m2)
 *
 *   pbm:   if m1 are overlaping or including into m2, this is formule will not be valide
 *   example AByyyCD, if sm = yyy, m1 = ByyyC, m2 = AByyyCD
 *   like that, sm exit in m1, we extract the nboccrs, after that we do the same with m2
 *   so, risk to have negative value, because we extrct the same value twice (the same positions)
 *
 *  solution:
 *  a) delete all the sub motif. and on the resulting matrix update.
 *  let the folowing example:
 *  s = AByyyCDzzzzByyyCzzzzyyy
 *  nb(yyy) = 3
 *  nb(ByyyC) = 2
 *  nb(AByyyCD) = 1
 *  on this example, (a) dosn't work, because all the motifs remain afeter deletion.
 *  and it will cause the pbm.
 *
 *  b)
 *  it may work, if we begin by the motif that have the longest size,
 *  and we search for its sm. and we decrease
 *  NO NO NO we fall in the same case and problem explain before
 *
 *  C) i could juste reduce by:
 *  look at the motif m1,m2, that have the greatest number of occrs, and substract that from the subMotifs
 *  this will better if we leave it whithout dcreasing.
 *
 *  D) i should look to the tree, in the moment of extracting the motif
 *  maybe i can find the solution, or maybe reduce the problem.
 */
void Old_CommonMotifs::delete_upDate_SubMotifsFromTheMatrix()
{
    cout<<" for now, this methode dosn't work, we need only the comments before it "<<endl;
}

unordered_map<uint32_t , uint32_t > Old_CommonMotifs::groupeCountSeqsByFamily(const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs)
{
    unordered_map<uint32_t , uint32_t > umap_familyId_ndSeqs{};

    uint32_t global_seqId;
    uint32_t family_id;

    for (const auto & pair_cm_seqId_nbOcc : cm_umap_seqId_nbOccs)
    {
        global_seqId = pair_cm_seqId_nbOcc.first;
        family_id = this->gst.get_FamilyId_And_SeqId(global_seqId).first; //(idx_family,idx_seq_in_family)

        if(umap_familyId_ndSeqs.find(family_id)==umap_familyId_ndSeqs.end())
        {
            umap_familyId_ndSeqs.insert(make_pair(family_id,1)); // the first seq in the family
        }
        else
        {
            umap_familyId_ndSeqs[family_id]++; // other seq in the same family
        }
    }

    return umap_familyId_ndSeqs;
}

bool Old_CommonMotifs::isCm_Common_For_ListSeqId_Percentage(
        const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs, uint32_t percentage_same_family,
        uint32_t percentage_multiple_families, uint32_t alpha)
{
    return false;
    /*
    to be continued ... check if same family ....
    uint32_t cm_nb_rep=0;

    uint32_t nb_all_seqs = 0;
    uint32_t nb_seqs_have_cm = 0;

    vector<int> list_nb_repeated_cm_by_seqId;

    for (const auto &seqId_nbOccs : cm_umap_seqId_nbOccs)
    {
        cm_nb_rep = seqId_nbOccs.second;
        list_nb_repeated_cm_by_seqId.push_back(cm_nb_rep);
    }

    pair<int,vector<int>> best_key_subSt = SubSetDistancePercentage::getBestSubSet(list_nb_repeated_cm_by_seqId,alpha);

    nb_seqs_have_cm = best_key_subSt.second.size();

    return ( (double)(nb_seqs_have_cm*100)/(double)nb_all_seqs >= percentage_same_family );

     */
}

void Old_CommonMotifs::generateMatrixCmsSeqs_deleteSubMotifs()
{
    string root_motif; // at the root, we have an empty motif.

    const Node * root = this->gst.getRootSuffixTree();

    for (auto & child : root->listChildren)
    {
        if(child)
        {
            fillMatrixCommonMotifsSeqsId_deleteSubMotifs(child, root_motif);
        }
    }
}

umap_seqId_nbOccs_ Old_CommonMotifs::fillMatrixCommonMotifsSeqsId_deleteSubMotifs(const Node *node, const string &prefixe)
{
    string motif_active_node = prefixe + node->cara;
    umap_seqId_nbOccs_ umap_seqId_nbOccs_active_node{};

    //if(motif_active_node.size() >= this->gst.getMinLengthMotif())
    if(motif_active_node.size() > this->gst.getMinLengthMotif())
    {
        if(!node->map_seqId_nbOccs.empty())
            umap_seqId_nbOccs_active_node = node->map_seqId_nbOccs; // COPY

        umap_seqId_nbOccs_ umap_seqId_nbOccs_sub_tree;

        for (auto & child : node->listChildren)
        {
            if(child)
            {
                umap_seqId_nbOccs_sub_tree = fillMatrixCommonMotifsSeqsId(child, motif_active_node);
                umapMerge(umap_seqId_nbOccs_active_node, umap_seqId_nbOccs_sub_tree);
            }
        }

        // at this stage, we have the umap_seqId_nbOccs_active_node that contain all seqId and position of the actuel node,
        // i.e: where the motifs that end in this node
        // add this common motifs to the matrix.

        //if(umap_seqId_nbOccs_active_node.size()>GAMMA)
        //if(!umap_seqId_nbOccs_active_node.empty() && isCm_Common_For_ListSeqId_Percentage(umap_seqId_nbOccs_active_node,BETA,BETA_2))
        //if(!umap_seqId_nbOccs_active_node.empty() && isCm_Common_For_ListSeqId_Percentage_cmInter(umap_seqId_nbOccs_active_node,BETA))
        if(!umap_seqId_nbOccs_active_node.empty() && isCm_Common_For_ListSeqId_Percentage_cmInter_2(umap_seqId_nbOccs_active_node,Beta))
        {
            // if the motif_active_node is a subMotif of already insrted motifs, we don't add it.
            // if it is not:
            // we add it to the matrix, and we check if therse motifs in the matrix
            // that are sub-motifs of this new insrted motif_active_node.
            if(!is_cm_sub_motifs(motif_active_node,umap_seqId_nbOccs_active_node))
            {
                // delete subMotifs of the motif_active_node
                delete_subMotifs(motif_active_node,umap_seqId_nbOccs_active_node);

                list_of_all_commonMotifs.push_back(motif_active_node);
                uint32_t idx_new_cm = list_of_all_commonMotifs.size()-1;
                map_str_len_idx.insert(make_pair(motif_active_node.length(),idx_new_cm));
                //bitVec_motifs.push_back(true);

                insert_cm_umap_seqId_nbOccs_In_matrix_commonMotifs_seqsId(umap_seqId_nbOccs_active_node);

                // important:
                // this code is important , not for here, but whene  i will construct the Htree directly without distMatrix, so i will need only the seqs that share the CM.
                //list_cm_umap_seqId_nbOccs.push_back(umap_seqId_nbOccs_active_node);
            }
        }
    }
    else
    {
        for (auto & child : node->listChildren)
        {
            if(child)
            {
                fillMatrixCommonMotifsSeqsId(child, motif_active_node);
            }
        }
    }

    return umap_seqId_nbOccs_active_node;
}

void Old_CommonMotifs::insert_cm_umap_seqId_nbOccs_In_matrix_commonMotifs_seqsId_deleteSubMotifs(
        const unordered_map<unsigned int, unsigned int> &umap_seqId_nbOccs)
{
    vector<uint32_t > list_seqId_nbOccs(this->gst.getNbAllSequences());

    for(const auto &pair_seqId_nbOccs:umap_seqId_nbOccs)
    {
        uint32_t i = pair_seqId_nbOccs.first;
        uint32_t val = pair_seqId_nbOccs.second;

        list_seqId_nbOccs.at(i) = val;
    }

    this->matrix_commonMotifs_seqsId.push_back(list_seqId_nbOccs);
}

/**
 *
 * @param new_cm  the new found commn motif
 * @param umap_new_cm unordred map of seqId, nb occrs
 * @return if new_cm is sub-motif of one motif stored before (list_cm, matrix).
 */
bool Old_CommonMotifs::is_cm_sub_motifs(const string& new_cm, const umap_seqId_nbOccs_& umap_new_cm)
{
    for (auto it= map_str_len_idx.upper_bound(new_cm.length()); it!=map_str_len_idx.end() ; ++it)
    {
        //cout<<" ( len: "<< it->first<<", idx: "<<it->second<<"), ";

        //----------------------------------------------------------
        //1 check if the sm is submotif of m
        if(list_of_all_commonMotifs.at(it->second).find(new_cm) ==std::string::npos )
            return false;

        // at this stage, the new_cm is submotif of m.
        //2) check if they share the same sequnces && number of occurences

        return areCommonMotifsInSameSeqs_sameNbOccrs(umap_new_cm, it->second);
    }

    return false;
}

/**
 *
 * @param new_cm string, new found motif
 * @param umap_new_cm u map of seqId,nbOcrrs for new_cm
 *
 * check the motifs already exist in the (list_motifs,matrix) that have len<new_cm
 * if they submotif of new_cm
 * if yes: delete them, by only deleting thier entry from the map of index,len "map_str_len_idx"
 * at the end (not here, the end of algorithm), we just copy the element by index that stil exist in "map_str_len_idx"
 *
 */
void Old_CommonMotifs::delete_subMotifs(const string &new_cm, const umap_seqId_nbOccs_ &umap_new_cm)
{
    auto it_begin = map_str_len_idx.begin();
    //1 check if the sm is submotif of m
    if(new_cm.find(list_of_all_commonMotifs.at(it_begin->second)) != std::string::npos )
    {
        //2) check if they share the same sequnces && number of occurences
        if(areCommonMotifsInSameSeqs_sameNbOccrs(umap_new_cm, it_begin->second))
        {
            // delete the index from map
            map_str_len_idx.erase(it_begin->first);
        }
    }

    for (auto it = map_str_len_idx.lower_bound(new_cm.length()); it!=map_str_len_idx.begin() ; --it)
    {
        if(new_cm.find(list_of_all_commonMotifs.at(it->second)) != std::string::npos )
        {
            //2) check if they share the same sequnces && number of occurences
            if(areCommonMotifsInSameSeqs_sameNbOccrs(umap_new_cm, it->second))
            {
                // delete the index from map
                map_str_len_idx.erase(it->first);
            }
        }
    }
}
