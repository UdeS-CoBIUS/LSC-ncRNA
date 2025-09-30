//
// Created by ibra on 5/5/2020.
//

#include <numeric>
#include <algorithm>
#include <fstream>
#include "CommonMotifs.h"
#include "hasher_Container.h"
#include "SubSetDistancePercentage.h"
#include "SequenceIdManager.h"

CommonMotifs::CommonMotifs(const SuffixTree_QuadraticTime &gst, unsigned int sameFamilyPercentageThreshold)
    : gst(gst), sameFamilyPercentageThreshold(sameFamilyPercentageThreshold)
{
    //nb_total_cm=0;
    //nb_deleted_cm_from_list_cm=0;
    //nb_deleted_cm_not_accepted=0;
    //nb_unique_child = 0;

    vec_flist.resize(gst.getMaxLengthMotif()+1,forward_list<uint32_t>()); //+1 to simplify the access, to avoid -1 (str.lenght vs str.lenght-1)
    
    this->nb_all_sequences = gst.getNbAllSequences();
}

CommonMotifs::CommonMotifs(const SuffixTree_QuadraticTime &gst, unsigned int sameFamilyPercentageThreshold,
                                                     int occurrenceVariationTolerance)
                : gst(gst),
                    sameFamilyPercentageThreshold(sameFamilyPercentageThreshold),
                    occurrenceVariationTolerance(occurrenceVariationTolerance)
{
    //nb_total_cm=0;
    //nb_deleted_cm_from_list_cm=0;
    //nb_deleted_cm_not_accepted=0;
    //nb_unique_child = 0;

    vec_flist.resize(gst.getMaxLengthMotif()+1,forward_list<uint32_t>()); //+1 to simplify the access, to avoid -1 (str.lenght vs str.lenght-1)

    this->nb_all_sequences = gst.getNbAllSequences();
}

CommonMotifs::CommonMotifs(const SuffixTree_QuadraticTime &gst, unsigned int sameFamilyPercentageThreshold,
                                                     int occurrenceVariationTolerance,
                                                     unsigned int nbOccrs_allowed)
                : gst(gst),
                    sameFamilyPercentageThreshold(sameFamilyPercentageThreshold),
                    occurrenceVariationTolerance(occurrenceVariationTolerance),
                    nbOccrs_allowed(nbOccrs_allowed)
{
    vec_flist.resize(gst.getMaxLengthMotif()+1,forward_list<uint32_t>()); //+1 to simplify the access, to avoid -1 (str.lenght vs str.lenght-1)

    this->nb_all_sequences = gst.getNbAllSequences();
    
    this->list_sum_nb_seqs = gst.getListSumNbSeqs();

}


void CommonMotifs::cmsExtractionSelection()
{
    string root_motif; // at the root, we have an empty motif.

    const Node * root = this->gst.getRootSuffixTree();

    for (auto & child : root->listChildren)
    {
        if(child)
        {
            cmsExtractionSelection_(child, root_motif);
        }
    }

    //cout<<"nb_unique_child : "<<nb_unique_child<<endl;
}

umap_seqId_nbOccs_ CommonMotifs::cmsExtractionSelection_(const Node *node, const string &prefix)
{
    string motif_active_node = prefix + node->cara;
    umap_seqId_nbOccs_ umap_seqId_nbOccs_active_node{};

    uint32_t nb_children_local = 0; // just for debugining

    if(motif_active_node.size() > this->gst.getMinLengthMotif())
    {
        if(!node->map_seqId_nbOccs.empty())
            umap_seqId_nbOccs_active_node = node->map_seqId_nbOccs; // COPY

        umap_seqId_nbOccs_ umap_seqId_nbOccs_sub_tree;

        for (auto & child : node->listChildren)
        {
            if(child)
            {
                umap_seqId_nbOccs_sub_tree = cmsExtractionSelection_(child, motif_active_node);
                umapMerge(umap_seqId_nbOccs_active_node, umap_seqId_nbOccs_sub_tree);

                if(!umap_seqId_nbOccs_sub_tree.empty())
                    nb_children_local++;
            }
        }

        //if(nb_children_local == 1) nb_unique_child++;

        // at this stage, we have the umap_seqId_nbOccs_active_node that contain all seqId and position of the actuel node,
        // i.e: where the motifs that end in this node
        // add this common motifs to the matrix.

        if(is_cm_accepted_according_To_selection_parameters(motif_active_node,umap_seqId_nbOccs_active_node))
        {
            list_cms.push_back(motif_active_node);
            list_cm_umap_seqId_nbOccs.push_back(umap_seqId_nbOccs_active_node);

            // the cm and the list have the same position (idx)
            //insert_cm_umap_seqId_nbOccs_In_matrix_commonMotifs_seqsId(umap_seqId_nbOccs_active_node);
        }
    }
    else
    {
        for (auto & child : node->listChildren)
        {
            if(child)
            {
                cmsExtractionSelection_(child, motif_active_node);
            }
        }
    }

    return umap_seqId_nbOccs_active_node;
}

/**
 * A- for each motif (node that represnte cm in gst):
 * -1- extrcat cm from gst
 * -2- check if it is acceptable according to selection paramaters
 * -3- check if it is not sub-string of existing cm
 * -4- add it to accpted cm
 * -5- seacrh for substring comparing to this new cm to delete theme.
 *
 * B-update infos, by kepping only non dleted informations
 * C-delete Unnecessary intermediate Memory Usage
 */
void CommonMotifs::cmsExtractionSelectionDeletionCSMs()
{
    string root_motif; // at the root, we have an empty motif.

    const Node * root = this->gst.getRootSuffixTree();

    for (auto & child : root->listChildren)
    {
        if(child)
        {
            cmsExtractionSelectionDeletionCSMs_(child, root_motif);
        }
    }

    update_valid_lists_cms_infos();
    cleanUnnecessaryMemoryUsage();

//    cout<<"end cmsExtractionSelectionDeletionCSMs() "<<endl;
//
//    cout<<"nb_total_cm : "<<nb_total_cm<<endl;
//    cout<<"nb_deleted_cm_not_accepted : "<<nb_deleted_cm_not_accepted<<endl;
//    cout<<"nb_deleted_cm_from_list_cm : "<<nb_deleted_cm_from_list_cm<<endl;
//    cout<<"nb_unique_child : "<<nb_unique_child<<endl;
}

/**
 *A- for each motif (node that represnte cm in gst):
 * -1- extrcat cm from gst
 * -2- check if it is acceptable according to selection paramaters
 * -3- check if it is not sub-string of existing cm
 * -4- add it to accpted cm
 * -5- seacrh for substring comparing to this new cm to delete theme.
 *
 * @param node
 * @param prefix
 * @return
 */
umap_seqId_nbOccs_ CommonMotifs::cmsExtractionSelectionDeletionCSMs_(const Node *node, const string &prefix)
{
    string motif_active_node = prefix + node->cara;
    umap_seqId_nbOccs_ umap_seqId_nbOccs_active_node{};

    //--------------------------------------------------------------------
    /*
    if(motif_active_node == "CCC")
    {
        cout<<"\n just for debugging porpose (at begin): ...."<<endl;
        //2) print for each paire: seqId and nbOccs:
        for (const auto & pair_seqId_nbOccs:umap_seqId_nbOccs_active_node)
        {
            cout<<"[ ";
            cout << pair_seqId_nbOccs.first; // seqId
            cout<<" : ";
            cout<<util::to_string(pair_seqId_nbOccs.second); // nbOccs
            cout<<" ]";
        }
    }
    */
    //--------------------------------------------------------------------

    if(motif_active_node.size() > this->gst.getMinLengthMotif())
    {
        if(!node->map_seqId_nbOccs.empty())
            umap_seqId_nbOccs_active_node = node->map_seqId_nbOccs; // COPY

        //--------------------------------------------------------------------
        /*
        if(motif_active_node == "CCC")
        {
            cout<<"\n just for debugging porpose (after if map is empty): ...."<<endl;
            //2) print for each paire: seqId and nbOccs:
            for (const auto & pair_seqId_nbOccs:umap_seqId_nbOccs_active_node)
            {
                cout<<"[ ";
                cout << pair_seqId_nbOccs.first; // seqId
                cout<<" : ";
                cout<<util::to_string(pair_seqId_nbOccs.second); // nbOccs
                cout<<" ]";
            }
        }
         */
        //--------------------------------------------------------------------

        umap_seqId_nbOccs_ umap_seqId_nbOccs_sub_tree;

        for (auto & child : node->listChildren)
        {
            if(child)
            {
                umap_seqId_nbOccs_sub_tree = cmsExtractionSelectionDeletionCSMs_(child, motif_active_node);

                //--------------------------------------------------------------------
                /*
                if(motif_active_node == "CCC")
                {
                    cout<<"\n just for debugging porpose (child subTree): ...."<<child->cara<<endl;
                    //2) print for each paire: seqId and nbOccs:
                    for (const auto & pair_seqId_nbOccs:umap_seqId_nbOccs_sub_tree)
                    {
                        cout<<"[ ";
                        cout << pair_seqId_nbOccs.first; // seqId
                        cout<<" : ";
                        cout<<util::to_string(pair_seqId_nbOccs.second); // nbOccs
                        cout<<" ]";
                    }
                }
                 */
                //--------------------------------------------------------------------

                umapMerge(umap_seqId_nbOccs_active_node, umap_seqId_nbOccs_sub_tree);

                //--------------------------------------------------------------------
                /*
                if(motif_active_node == "CCC")
                {
                    cout<<"\n just for debugging porpose (after Merge child): ...."<<child->cara<<endl;
                    //2) print for each paire: seqId and nbOccs:
                    for (const auto & pair_seqId_nbOccs:umap_seqId_nbOccs_active_node)
                    {
                        cout<<"[ ";
                        cout << pair_seqId_nbOccs.first; // seqId
                        cout<<" : ";
                        cout<<util::to_string(pair_seqId_nbOccs.second); // nbOccs
                        cout<<" ]";
                    }
                }
                 */
                //--------------------------------------------------------------------
            }
        }

        //--------------------------------------------------------------------
        /*
        if(motif_active_node == "CCC")
        {
            cout<<"\n just for debugging porpose (after sub-tree): ...."<<endl;
            //2) print for each paire: seqId and nbOccs:
            for (const auto & pair_seqId_nbOccs:umap_seqId_nbOccs_active_node)
            {
                cout<<"[ ";
                cout << pair_seqId_nbOccs.first; // seqId
                cout<<" : ";
                cout<<util::to_string(pair_seqId_nbOccs.second); // nbOccs
                cout<<" ]";
            }
        }
         */
        //--------------------------------------------------------------------

        // at this stage, we have the umap_seqId_nbOccs_active_node that contain all seqId and position of the actuel node,
        // i.e: where the motifs that end in this node
        // add this common motifs to the matrix.

        if(is_cm_accepted_according_To_selection_parameters(motif_active_node,umap_seqId_nbOccs_active_node))
        {
            //nb_total_cm++;

            std::map<unsigned int , unsigned int > temp_map(umap_seqId_nbOccs_active_node.begin(), umap_seqId_nbOccs_active_node.end());
            auto ncm_hashValue = hasher_Container::hash_64bits(temp_map);

            if(!is_ncm_sub_motif_list_cms(motif_active_node,umap_seqId_nbOccs_active_node,ncm_hashValue))
            {

                // -------------------------------------------------------------
                /*
                if(motif_active_node == "CCC")
                {
                    cout<<" just for debugging porpose (adding to the matrix): ...."<<endl;
                    //2) print for each paire: seqId and nbOccs:
                    for (const auto & pair_seqId_nbOccs:umap_seqId_nbOccs_active_node)
                    {
                        cout<<"[ ";
                        cout << pair_seqId_nbOccs.first; // seqId
                        cout<<" : ";
                        cout<<util::to_string(pair_seqId_nbOccs.second); // nbOccs
                        cout<<" ]";
                    }
                }
                */
                // -------------------------------------------------------------

                add_ncm_to_lists_cms(motif_active_node,umap_seqId_nbOccs_active_node,ncm_hashValue);

                delete_sub_motifs_of_ncm_from_list_cms(motif_active_node, umap_seqId_nbOccs_active_node, ncm_hashValue);
            }
            else{//nb_deleted_cm_not_accepted++;
                }//this for debuging only
        }
    }
    else
    {
        for (auto & child : node->listChildren)
        {
            if(child)
            {
                cmsExtractionSelectionDeletionCSMs_(child, motif_active_node);
            }
        }
    }

    return umap_seqId_nbOccs_active_node;
}


/**
 * delete the submotif from the matrix and the list of motifs as well
 * we delete if: 1) the sm and m share the same sequences and
 *               2) have the same number of occurences
 */
void CommonMotifs::deleteSubMotifs_at_end_ES()
{
    unsigned int nb_deleted_elements=0;

    // 1) creat an index of all motifs, and sort it according to the motis lenght
    vector<int > list_motif_id(list_cms.size());
    std::iota(list_motif_id.begin(),list_motif_id.end(),0);

    std::sort(list_motif_id.begin(),list_motif_id.end(),[&](uint32_t i, uint32_t j){
        return list_cms.at(i).size() < list_cms.at(j).size();
    });

    // only check for debug
    //cout<<"before deleting : "<<endl;
    //for(unsigned int i =0; i<list_motif_id.size(); ++i) cout<<"| "<<list_of_all_commonMotifs.at(list_motif_id.at(i));

    // 2) check if motis i is sub-motifs of one of motifs that are greater than it. if yes delete it.
    for(unsigned int i =0; i<list_motif_id.size()-1;++i)
    {
        for (unsigned int j = i+1; j < list_motif_id.size(); ++j)
        {

            if(is_sub_motif(list_motif_id.at(i), list_motif_id.at(j)))
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
    vector<string> tmp_list_cms;
    tmp_list_cms.reserve(list_cms.size()-nb_deleted_elements);

    vector< umap_seqId_nbOccs_ > tmp_list_cm_umap_seqId_nbOccs;
    tmp_list_cm_umap_seqId_nbOccs.reserve(list_cm_umap_seqId_nbOccs.size()-nb_deleted_elements);

    for(int i : list_motif_id){
        if(i!=-1){
            tmp_list_cms.push_back( list_cms.at(i) );
            tmp_list_cm_umap_seqId_nbOccs.push_back( list_cm_umap_seqId_nbOccs.at(i) );
        }
    }

    //cout<<"\n tmp-list_strings : "<<endl;
    //cout<<util::to_string(tmp_list_of_all_commonMotifs);

    this->list_cms = tmp_list_cms;
    this->list_cm_umap_seqId_nbOccs = tmp_list_cm_umap_seqId_nbOccs;

    cout<<"nb motifs after deletion : "<<list_cms.size()<<endl;
}

/**
 * delete the submotif from the matrix and the list of motifs as well
 * we delete if: 1) the sm and m share the same sequences and
 *               2) have the same number of occurences
 *
 * using hashage to accelerating the comparaison if each column
 * we use tow hash table
 * 1) store the hash value only for seqId
 * 2) store the hash value for nb_occrs
 *
 * so in the test, we check 1) if the hashvalue is the same for seqsId,
 * then we go to check if the hash vlue is the same for nb_occrs
 *
 * ofcours, the best is to check only once, by creating only one hash table
 * that combine the 1 and 2. it's depend on the hash function be used.
 *
 */
void CommonMotifs::deleteSubMotifs_at_end_ES_usingHash()
{
    unsigned int nb_deleted_elements=0;

    // 1) creat an index of all motifs, and sort it according to the motifs length
    vector<int > list_motif_id(list_cms.size());
    std::iota(list_motif_id.begin(),list_motif_id.end(),0);

    std::sort(list_motif_id.begin(),list_motif_id.end(),[&](uint32_t i, uint32_t j){
        return list_cms.at(i).size() < list_cms.at(j).size();
    });

    // 1") creat hash value for the umap
    list_hashValue_cm_umap.reserve(list_cm_umap_seqId_nbOccs.size());

    for (auto & cm_umap : list_cm_umap_seqId_nbOccs)
    {
        std::map<uint32_t ,uint32_t > temp_map(cm_umap.begin(),cm_umap.end());
        list_hashValue_cm_umap.push_back(hasher_Container::hash_64bits(temp_map));
    }


    // only check for debug
    //cout<<"before deleting : "<<endl;
    //for(unsigned int i =0; i<list_motif_id.size(); ++i) cout<<"| "<<list_of_all_commonMotifs.at(list_motif_id.at(i));

    // 2) check if motis i is sub-motifs of one of motifs that are greater than it. if yes delete it.
    for(unsigned int i =0; i<list_motif_id.size()-1;++i)
    {
        for (unsigned int j = i+1; j < list_motif_id.size(); ++j)
        {
            //bool can_delete = is_sub_motif(list_cms.at(i),list_cms.at(j));
            //cout<<""

                //1 check if the sm is submotif of lm
            //2) check if they share the same sequnces && the number of occurences
            if(
                    list_cms.at(list_motif_id.at(j)).find(list_cms.at(list_motif_id.at(i))) !=std::string::npos
                    &&
                    list_hashValue_cm_umap.at(list_motif_id.at(i)) == list_hashValue_cm_umap.at(list_motif_id.at(j))
                    //&& // this last condition if hashValue are equale, we check to avoid colision
                    //list_cm_umap_seqId_nbOccs.at(list_motif_id.at(i)) == list_cm_umap_seqId_nbOccs.at(list_motif_id.at(j))
                    )
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
    vector<string> tmp_list_cms;
    tmp_list_cms.reserve(list_cms.size()-nb_deleted_elements);

    vector< umap_seqId_nbOccs_ > tmp_list_cm_umap_seqId_nbOccs;
    tmp_list_cm_umap_seqId_nbOccs.reserve(list_cm_umap_seqId_nbOccs.size()-nb_deleted_elements);

    for(int i : list_motif_id){
        if(i!=-1){
            tmp_list_cms.push_back( list_cms.at(i) );
            tmp_list_cm_umap_seqId_nbOccs.push_back( list_cm_umap_seqId_nbOccs.at(i) );
        }
    }

    //cout<<"\n tmp-list_strings : "<<endl;
    //cout<<util::to_string(tmp_list_of_all_commonMotifs);

    this->list_cms = tmp_list_cms;
    this->list_cm_umap_seqId_nbOccs = tmp_list_cm_umap_seqId_nbOccs;

    cout<<"nb motifs after deletion : "<<list_cms.size()<<endl;
}


/**
 * check if sm is real sub motif of lm, so we can delete it
 * 1) if sm is a sub-string of lm
 * 2) sm && lm have the same list of seqId and nbOccrs :
 *   they exist in the same sequences and have the same nbOcrrs in each sequence.
 * @param idx_sm idx of small motif in list_cms && list_cm_umap_seqId_nbOccs  |sm|<=|lm|
 * @param idx_m idx of longer motif in list_cms && list_cm_umap_seqId_nbOccs |sm|<=|lm|
 * @return if we can delete the submotifs
 */
inline bool CommonMotifs::is_sub_motif(unsigned int idx_sm, unsigned int idx_m)
{
    //1 check if the sm is submotif of lm
    //2) check if they share the same sequnces && the number of occurences

    return (
           list_cms.at(idx_m).find(list_cms.at(idx_sm)) != std::string::npos
           &&
           list_cm_umap_seqId_nbOccs.at(idx_sm) == list_cm_umap_seqId_nbOccs.at(idx_m)
            );
}

inline bool CommonMotifs::is_sub_motif_hash(unsigned int idx_sm, unsigned int idx_m)
{
    return (
            list_cms.at(idx_m).find(list_cms.at(idx_sm)) != std::string::npos
            &&
            list_hashValue_cm_umap.at(idx_sm) == list_hashValue_cm_umap.at(idx_m)
            //&&
            //list_cm_umap_seqId_nbOccs.at(idx_sm) == list_cm_umap_seqId_nbOccs.at(idx_m)
    );
}

bool CommonMotifs::is_sub_motif_hash(const string &ncm, const umap_seqId_nbOccs_ &ncm_umap, uint64_t ncm_hashValue, uint32_t idx_m)
{
    return(
            list_cms.at(idx_m).find(ncm) != std::string::npos
            &&
            list_hashValue_cm_umap.at(idx_m) == ncm_hashValue
            &&
            ncm_umap==list_cm_umap_seqId_nbOccs.at(idx_m)
            );
}

bool CommonMotifs::is_ncm_super_motif_hash(const string &ncm, const umap_seqId_nbOccs_ &ncm_umap,
                                           uint64_t ncm_hashValue, uint32_t idx_m)
{
    return(
            ncm.find(list_cms.at(idx_m)) != std::string::npos
            &&
            list_hashValue_cm_umap.at(idx_m) == ncm_hashValue
            //&&
            //ncm_umap==list_cm_umap_seqId_nbOccs.at(idx_m)
    );
}

/**
 * get list_motifs length > new_motif
 * check if it sub-motif, return true
 * at the end of the list return true
 *
 * @param ncm  the new found commn motif
 * @param ncm_umap unordred map of seqId, nb occrs
 * @return if new_cm is sub-motif of one motif stored before (list_cm, matrix).
 */
bool
CommonMotifs::is_ncm_sub_motif_list_cms(const string &ncm, const umap_seqId_nbOccs_ &ncm_umap, uint64_t ncm_hashValue)
{

    uint32_t i = ncm.length()+1; // upper bound, to check with string that have > length
    for(;i<vec_flist.size();++i)
    {
        for (uint32_t idx : vec_flist.at(i))
        {
            if(is_sub_motif_hash(ncm, ncm_umap, ncm_hashValue, idx)) return true;
        }
    }

    return false;
}

void
CommonMotifs::add_ncm_to_lists_cms(string &ncm, umap_seqId_nbOccs_ &ncm_umap, unsigned long long int ncm_hashValue)
{
    list_cms.push_back(ncm);
    list_cm_umap_seqId_nbOccs.push_back(ncm_umap);
    list_hashValue_cm_umap.push_back(ncm_hashValue);

    // add ncm idx in list_cm to the vec_flist
    vec_flist.at(ncm.length()).push_front(list_cms.size()-1);
}

/**
 * if we find cm in list_cm that are subMotifs of ncm, we delete theme.
 * @param ncm
 * @param ncm_umap
 * @param ncm_hashValue
 * @param flist
 */
void CommonMotifs::flist_erase_all_idx_of_subMotif_of_ncm(string &ncm, umap_seqId_nbOccs_ &ncm_umap,
                                                          unsigned long long int ncm_hashValue, forward_list<uint32_t> &flist)
{

    auto prev = flist.before_begin();
    for (auto it = flist.begin(); it!=flist.end(); )
    {
        if(is_ncm_super_motif_hash(ncm, ncm_umap, ncm_hashValue, *it))
        {
            it = flist.erase_after(prev);
            // break; // comment or uncomment to deal with only the first or all found element(s).

            //nb_deleted_cm_from_list_cm++;
        }
        else
        {
            prev = it;
            ++it;
        }
    }
}

/**
 * ---- get list_motifs length < new_motif
---- for each motif_in_list check if its sub-motif of new_motif
---- yes delete it: from the map of valid_motif

 * @param ncm
 * @param ncm_umap
 * @param ncm_hashValue
 */
void CommonMotifs::delete_sub_motifs_of_ncm_from_list_cms(string &ncm, umap_seqId_nbOccs_ &ncm_umap,
                                                          unsigned long long int ncm_hashValue)
{
    // check all motifs that are < ncm
    // we begin from MinLengthMotif, because before there is no motifs
    // to : ncm.length() -1
    for(uint32_t i = gst.getMinLengthMotif(); i < ncm.length(); ++i)
    {
        // for each flist in vec_flist : delete all submotifs of ncm
        flist_erase_all_idx_of_subMotif_of_ncm(ncm, ncm_umap, ncm_hashValue, vec_flist.at(i));

    }
}

uint32_t CommonMotifs::get_nb_elements_in_vec_flist()
{
    uint32_t nb=0;

    for (auto &flist : vec_flist)
    {
        nb += std::distance(flist.begin(),flist.end());
    }

    return nb;
}

/**
 * update the matrix : keep only motifs that are not deleted in the matrix
 */
void CommonMotifs::update_valid_lists_cms_infos()
{
    //cout<<"begin update_valid_lists_cms_infos() : "<<map_str_len_idx.size()<<endl;

    uint32_t nb_remaining_elements = get_nb_elements_in_vec_flist();

    vector<string> tmp_list_cms;
    tmp_list_cms.reserve(nb_remaining_elements);

    vector< umap_seqId_nbOccs_ > tmp_list_cm_umap;
    tmp_list_cm_umap.reserve(nb_remaining_elements);


    for (auto &flist : vec_flist)
    {
        for (uint32_t i:flist)
        {
            tmp_list_cms.push_back( list_cms.at(i) );
            tmp_list_cm_umap.push_back( list_cm_umap_seqId_nbOccs.at(i) );
        }
    }

    this->list_cms = tmp_list_cms;
    this->list_cm_umap_seqId_nbOccs = tmp_list_cm_umap;
}

void CommonMotifs::cleanUnnecessaryMemoryUsage()
{
    //this->map_str_len_idx.clear();
    //this->list_hashValue_cm_umap.clear(); this->list_hashValue_cm_umap.shrink_to_fit();

    vector<forward_list<uint32_t >>().swap(this->vec_flist);
    vector<uint64_t >().swap(this->list_hashValue_cm_umap);
}


void CommonMotifs::generateMatrixCmsSeqs()
{
    uint32_t nb_total_seqs = this->nb_all_sequences;
    uint32_t nb_total_cms = this->list_cms.size(); // the same of : list_cm_umap_seqId_nbOccs.size();

    /*
    this->matrix_nbOcrrs_cmsSeqsIds.reserve(nb_total_cms);

    for (auto &cm_umap_seqId_nbOcc : this->list_cm_umap_seqId_nbOccs)
    {
        vector<uint32_t > list_seqId_nbOccs(nb_total_seqs);

        for(const auto &pair_seqId_nbOcc:cm_umap_seqId_nbOcc)
        {
            uint32_t i = pair_seqId_nbOcc.first;
            uint32_t val = pair_seqId_nbOcc.second;

            list_seqId_nbOccs.at(i) = val;
        }

        this->matrix_nbOcrrs_cmsSeqsIds.push_back(list_seqId_nbOccs);
    }
     */

    // resize the vector to R:nb_total_seqs elements of type vector<unsigned int>, each having size C:nb_total_cms
    this->matrix_nbOcrrs_cmsSeqsIds.resize(nb_total_seqs, std::vector<unsigned int>(nb_total_cms,0));

    for (unsigned int idx_cm=0; idx_cm < this->list_cm_umap_seqId_nbOccs.size(); ++idx_cm)
    {
        for(const auto &pair_seqId_nbOcc:list_cm_umap_seqId_nbOccs.at(idx_cm))
        {
            uint32_t idx_seq = pair_seqId_nbOcc.first;
            uint32_t nb_occrs = pair_seqId_nbOcc.second;

            this->matrix_nbOcrrs_cmsSeqsIds[idx_seq][idx_cm]=nb_occrs;
        }
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// ------------------------------------ CMs selection and filtrate Methods

bool CommonMotifs::is_cm_accepted_according_To_selection_parameters(const string &cm,
                                                                    const umap_seqId_nbOccs_ &cm_umap_seqId_nbOccs)
{
    //if(umap_seqId_nbOccs_active_node.size()>GAMMA)
    //if(!umap_seqId_nbOccs_active_node.empty() && is_cm_accepted_only_for_one_family_percentage_nbOcrrs(
    //        umap_seqId_nbOccs_active_node, sameFamilyPercentageThreshold, sameFamilyPercentageThresholdInternal))

    //return (!cm_umap_seqId_nbOccs.empty() && is_cm_accepted_only_for_one_family_percentage_nbOcrrs(cm_umap_seqId_nbOccs, sameFamilyPercentageThreshold));

    if(occurrenceVariationTolerance == -1) // only same-family percentage threshold
        return (cm_umap_seqId_nbOccs.size() >= GAMMA &&
                is_cm_accepted_only_for_one_family_percentage_nbOcrrs(cm_umap_seqId_nbOccs,
                                                                      sameFamilyPercentageThreshold,
                                                                      this->nbOccrs_allowed));

    // Otherwise check both thresholds plus nbOccrs_allowed
    return is_cm_accepted_only_for_one_family_percentage_and_variation(cm_umap_seqId_nbOccs,
                                                                       sameFamilyPercentageThreshold,
                                                                       occurrenceVariationTolerance,
                                                                       this->nbOccrs_allowed);
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
bool CommonMotifs::is_cm_accepted_only_for_one_family_percentage_nbOcrrs(
        const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs,
        uint32_t percentage_same_family,
        uint32_t nbOccrs_allowed_local)
{
    uint32_t nb_all_seqs_comparedTo = 0;
    uint32_t nb_seqs_have_cm = cm_umap_seqId_nbOccs.size();  // nb seqs that have the common motif.

    // if we use nbOccrs_allowed_local (which is > 1), we count all seqs that have the nbOccrs of the cm above or equqle the nbOccrs_allowed_local
    // todo: this should be donne at the time of adding nbOccrs to the maps, we check if it is above or equal to nbOccrs_allowed_local
    //       unless, we want to do the multiple parameter combination with the percentage threshold and nbOccrs like: nbOccrs=1 threshold = 40 and nbOccrs=2 threshold = 10 for example.
    if(nbOccrs_allowed_local>1) // 1 is default
    {
        uint32_t nb_seqs_have_nbOccrs_allowd = 0;

        for (auto pair_seqId_nbOccrs: cm_umap_seqId_nbOccs)
        {
            if (pair_seqId_nbOccrs.second >= nbOccrs_allowed_local)
                nb_seqs_have_nbOccrs_allowd++;
        }

        nb_seqs_have_cm = nb_seqs_have_nbOccrs_allowd;

        // if sameFamilyPercentageThreshold (percentage_same_family) = 0 and nb_seqs_have_cm<2,
        // this means that we have only one single seq that have the motif with nbOccrs_allowed
        // hence we don't accept it
        //if(nb_seqs_have_cm<2) // it dosen't metter if percentage_same_family==0 or its >0, we don't accepte
        if(percentage_same_family == 0 && nb_seqs_have_cm<2) // this mean for samll family a pbm, if nbsqs = 10, this mean mean threshold is = 20, if nbsqs = 4 threshold = 50, if nbsqs=3 threshold 70, nsqs=2 threshold = 100
            return false;
    }

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
    nb_all_seqs_comparedTo = this->nb_all_sequences;
    return ( (double)(nb_seqs_have_cm*100)/(double)nb_all_seqs_comparedTo >= percentage_same_family );
}

bool CommonMotifs::is_cm_accepted_only_for_one_family_percentage_and_variation(
    const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs,
    uint32_t percentage_same_family,
    uint32_t variation_tolerance,
    uint32_t nbOccrs_allowed_local)
{
    uint32_t nb_all_seqs_comparedTo = 0;

    // 1st: first, map the element of number of occurece of each seq in one single vector.
    uint32_t cm_nb_occrs=0;
    uint32_t nb_seqs_have_cm_with_alpha = 0; // nb of elements in the best list

    vector<int> list_nb_repeated_cm_by_seqId;

    for (const auto &seqId_nbOccs : cm_umap_seqId_nbOccs)
    {
        cm_nb_occrs = seqId_nbOccs.second;
        list_nb_repeated_cm_by_seqId.push_back(cm_nb_occrs);
    }

    // 2nd:  get the best list that satisfy the occurrence variation parameter and the same-family percentage threshold
    //      +
    //      the nbOccrs_allowed_local paramter if it is >= 2 , because, otherwise we will have 1 which is the defaut or the lowest value.

    pair<int,vector<int>> best_key_subSt;
    if(nbOccrs_allowed_local < 2)
        best_key_subSt = SubSetDistancePercentage::getBestSubSet(list_nb_repeated_cm_by_seqId, variation_tolerance);
    else
        best_key_subSt = SubSetDistancePercentage::getBestSubSet_lowerBound(list_nb_repeated_cm_by_seqId,
                                                                            variation_tolerance,
                                                                            nbOccrs_allowed_local);

    nb_seqs_have_cm_with_alpha = best_key_subSt.second.size();

    // 3rd: check if the cm satisafy the same-family percentage threshold; when it is 0 the motif must exist in at least 2 sequences.
    if(percentage_same_family == 0 && nb_seqs_have_cm_with_alpha >= GAMMA) {
    sum_cm_nb_occrs_within_variation += best_key_subSt.first;
        return true;
    } else if(percentage_same_family == 0 && nb_seqs_have_cm_with_alpha < GAMMA){
        return false;
    }

    if(this->gst.isSequencesAreGroupedByFamilies())
    {
        const auto umap_family_nbSeqs = groupeCountSeqsByFamily(cm_umap_seqId_nbOccs);

        for(const auto & pair_familyId_nbSeqs:umap_family_nbSeqs)
        {
            uint32_t first_family_id = pair_familyId_nbSeqs.first;
            nb_all_seqs_comparedTo = this->gst.getListFamiliesSeqs().at(first_family_id).size();

            if( ((double)(nb_seqs_have_cm_with_alpha*100)/(double)nb_all_seqs_comparedTo) >= percentage_same_family )
            {
                sum_cm_nb_occrs_within_variation += best_key_subSt.first;
                return true;
            }
        }

        return false;
    }

    // here, meaning sequnces are not grouped, so they are supposed tu be all in the same family
    nb_all_seqs_comparedTo = this->nb_all_sequences;

    if ( (double)(nb_seqs_have_cm_with_alpha*100)/(double)nb_all_seqs_comparedTo >= percentage_same_family )
    {
    sum_cm_nb_occrs_within_variation += best_key_subSt.first;
        return true;
    }

    return false;
}

// ---------------------------------------------------------------------------------------------------------------------
// ------------------------------------ Utils Methods

void CommonMotifs::umapMerge(umap_seqId_nbOccs_ &des, umap_seqId_nbOccs_ &src)
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

unordered_map<uint32_t , uint32_t > CommonMotifs::groupeCountSeqsByFamily(const unordered_map<unsigned int, unsigned int> &cm_umap_seqId_nbOccs)
{
    unordered_map<uint32_t , uint32_t > umap_familyId_ndSeqs{};

    uint32_t global_seqId;
    uint32_t family_id;

    for (const auto & pair_cm_seqId_nbOcc : cm_umap_seqId_nbOccs)
    {
        global_seqId = pair_cm_seqId_nbOcc.first;
        family_id = this->get_FamilyId_And_SeqId(global_seqId).first; //(idx_family,idx_seq_in_family)

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

void CommonMotifs::print_list_cm_umap_seqId_nbOccs() const
{
    for (uint32_t i=0 ; i<list_cm_umap_seqId_nbOccs.size(); ++i)
    {
        const auto & motif = list_cms.at(i);
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

void CommonMotifs::print_matrix_cms_seqsId() const
{
    unsigned int nb_seqs = nb_all_sequences;
    unsigned int nb_motifs = this->list_cms.size();

    // 1) print the motifs
    cout<<endl<<"        ";
    for (const auto & motif:this->list_cms)
    {
        cout<<motif <<"|  ";
    }
    cout<<endl;

    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        cout<<"seq("<<i<<"): ";

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            cout<<this->matrix_nbOcrrs_cmsSeqsIds[i][j];
            cout<<"| ";
        }

        cout<<"\n____________"<<endl;
    }
}

void CommonMotifs::print_infos() const
{
    std::cout << "Printing CommonMotifs infos:" << std::endl;

    unsigned int nb_seqs = this->nb_all_sequences;
    unsigned int nb_motifs = this->list_cms.size();

    cout<<" Total nb_seqs = "<<nb_seqs<<endl;
    cout<<" Total nb_motifs = "<<nb_motifs<<endl;
    cout<<" sum_cm_nb_occrs_within_variation = "<<sum_cm_nb_occrs_within_variation<<endl;
    cout<<" Average nb_occres = "<<(double (sum_cm_nb_occrs_within_variation) / double (nb_motifs))<<endl;
    std::cout << "Number of rows (sequences): " << matrix_nbOcrrs_cmsSeqsIds.size() << std::endl;
    std::cout << "Number of columns (common motifs): " << (matrix_nbOcrrs_cmsSeqsIds.empty() ? 0 : matrix_nbOcrrs_cmsSeqsIds[0].size()) << std::endl;
}

// using buffer , and then use file "myfilestream << buffer".
void CommonMotifs::saveMatrixCMS_ToCsv_File(string file_output) const
{
    unsigned int nb_seqs= this->nb_all_sequences;
    unsigned int nb_motifs = this->list_cms.size();

    // nb seqs in the matrix_nbOcrrs_cmsSeqsIds
    nb_seqs = matrix_nbOcrrs_cmsSeqsIds.size(); // this work fine, the other I don't know why it give 1, and in other function it work fine
   
    ofstream file_csv;
    const string csv_separator =",";
    string buffer;
    uint32_t size_buffer = 0;

    size_buffer += nb_motifs * nb_seqs; // the size of the matrix_commonMotifs_seqsId

    // this supose that we have only one digit, but we can have tow digit!! so howe  i should do this
    // generaly we have only one digit, but tow digit exist.
    // i think to add 1/3 of the space to cover the possible tow digit
    size_buffer += size_buffer/3;


    size_buffer += 2*nb_seqs; // we add 1 colomn for familyId, we assum that we have 2 degit by each familyId
    size_buffer += nb_seqs; // each seq is a row, so it is an "\n" at the end

    // 1) familyId is the class.

    uint32_t size_all_colomn_name = string("familyId").size();


    // 2) compute the total size of all motifs
    for (const auto & motif:this->list_cms)
    {
        size_all_colomn_name+=motif.size()+1; // +1 for the separtor
    }
    size_all_colomn_name +=1; // for the endline of column name

    size_buffer+=size_all_colomn_name;

    // reserve the space memory , then add the informations
    buffer.reserve(size_buffer);

    // the column names:
    buffer+="familyId";
    for (const auto & motif:this->list_cms)
    {
        buffer+=csv_separator;
        buffer+=motif;
    }
    buffer+="\n";

    // the matrix values
    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->gst.isSequencesAreGroupedByFamilies())
        {
            //std::cout << "this->gst.isSequencesAreGroupedByFamilies() = True, i: " << i << std::endl;
            auto pair_familyId_SeqId = this->get_FamilyId_And_SeqId(i); //(idx_family,idx_seq_in_family)
            //std::cout << "pair_familyId_SeqId.first: " << pair_familyId_SeqId.first << std::endl;
            buffer+=util::to_string(pair_familyId_SeqId.first); //idx_family
        }
        else
        {
            buffer+="0"; //idx_family is 0, we have one family
        }

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            buffer+=csv_separator;
            buffer+=util::to_string(matrix_nbOcrrs_cmsSeqsIds[i][j]); // nb occurence
        }
        buffer+="\n";
    }

    //-------------------------------------------------------------------
    // string dir_path = file_output+".csv";
    string dir_path = file_output; // the file_output already have the .csv extension
    file_csv.open (dir_path);
    file_csv<<buffer;

    file_csv.close();

    cout << "Matrix has been saved to " << dir_path << endl;
}


// using buffer , and then use file "myfilestream << buffer".
// use std::to_string instead of util::to_string (my implemenation use ostringstream)
void CommonMotifs::saveMatrixCMS_ToCsv_File_stdtostring(string file_output) const
{
    unsigned int nb_seqs= this->nb_all_sequences;
    unsigned int nb_motifs = this->list_cms.size();

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


    size_buffer += 2*nb_seqs; // we add 1 colomn for familyId, we assum that we have 2 degit by each familyId
    size_buffer += nb_seqs; // each seq is a row, so it is an "\n" at the end

    // 1) familyId is the class.

    uint32_t size_all_colomn_name = string("familyId").size();


    // 2) compute the total size of all motifs
    for (const auto & motif:this->list_cms)
    {
        size_all_colomn_name+=motif.size()+1; // +1 for the separtor
    }
    size_all_colomn_name +=1; // for the endline of column name

    size_buffer+=size_all_colomn_name;

    // reserve the space memory , then add the informations
    buffer.reserve(size_buffer);

    // the column names:
    buffer+="familyId";
    for (const auto & motif:this->list_cms)
    {
        buffer+=csv_separator;
        buffer+=motif;
    }
    buffer+="\n";

    // the matrix values
    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->gst.isSequencesAreGroupedByFamilies())
        {
            auto pair_familyId_SeqId = this->get_FamilyId_And_SeqId(i); //(idx_family,idx_seq_in_family)
            buffer+=std::to_string(pair_familyId_SeqId.first); //idx_family
        }
        else
        {
            buffer+="0"; //idx_family is 0, we have one family
        }

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            buffer+=csv_separator;
            buffer+=std::to_string(matrix_nbOcrrs_cmsSeqsIds[i][j]); // nb occurence
        }
        buffer+="\n";
    }

    //-------------------------------------------------------------------
    string dir_path = file_output+".csv";
    file_csv.open (dir_path);
    cout<<" file_csv.open (\"matrix_motifs_seqIds.csv\"); "<<endl;

    file_csv<<buffer;

    file_csv.close();
    cout<<"file_csv.close();"<<endl;
}


// using buffer , and then use file "myfilestream << buffer".
// use std::to_string instead of util::to_string (my implemenation use ostringstream)
// don't compute the total size by going through each element.
// use the chatGPT idea form the other function, just give an istimated size
void CommonMotifs::saveMatrixCMS_ToCsv_File_stdtostring_size_estimated(string file_output) const
{
    unsigned int nb_seqs= this->nb_all_sequences;
    unsigned int nb_motifs = this->list_cms.size();

    cout<<" Total nb_seqs = "<<nb_seqs<<endl;
    cout<<" Total nb_motifs = "<<nb_motifs<<endl;

    ofstream file_csv;
    const string csv_separator =",";
    string buffer;
    uint32_t size_buffer = 0;


    // the header information size:
    size_buffer = nb_motifs*10 + nb_motifs + 20; // we suppse each motif take in average 10 (mid [2-20]) + nb_motifs for separtor, + 20 (suplimenatry and include siezof('familyId'))

    // the matrix
    size_buffer += nb_motifs * nb_seqs * 4; // (nb_motifs * nb_seqs) the size of the matrix_commonMotifs_seqsId, *4 supose each digit take 4
    size_buffer += nb_seqs*4; // we add 1 colomn for familyId = nb_seqs entry, *4 : each digit take 4
    size_buffer += nb_seqs; // each seq is a row, so it is an "\n" at the end
    
    // reserve the space memory , then add the informations
    buffer.reserve(size_buffer);

    // the column names:
    buffer+="familyId";
    for (const auto & motif:this->list_cms)
    {
        buffer+=csv_separator;
        buffer+=motif;
    }
    buffer+="\n";

    // the matrix values
    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->gst.isSequencesAreGroupedByFamilies())
        {
            auto pair_familyId_SeqId = this->get_FamilyId_And_SeqId(i); //(idx_family,idx_seq_in_family)
            buffer+=std::to_string(pair_familyId_SeqId.first); //idx_family
        }
        else
        {
            buffer+="0"; //idx_family is 0, we have one family
        }

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            buffer+=csv_separator;
            buffer+=std::to_string(matrix_nbOcrrs_cmsSeqsIds[i][j]); // nb occurence
        }
        buffer+="\n";
    }

    //-------------------------------------------------------------------
    string dir_path = file_output+".csv";
    file_csv.open (dir_path);
    cout<<" file_csv.open (\"matrix_motifs_seqIds.csv\"); "<<endl;

    file_csv<<buffer;

    file_csv.close();
    cout<<"file_csv.close();"<<endl;
}

// chatGPT: use  `std::ofstream::write` function instead of `operator<<`
void CommonMotifs::saveMatrixCMS_ToCsv_File_ofstream_write(string file_output) const
{
    unsigned int nb_seqs= this->nb_all_sequences;
    unsigned int nb_motifs = this->list_cms.size();

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


    size_buffer += 2*nb_seqs; // we add 1 colomn for familyId, we assum that we have 2 degit by each familyId
    size_buffer += nb_seqs; // each seq is a row, so it is an "\n" at the end

    // 1) familyId is the class.

    uint32_t size_all_colomn_name = string("familyId").size();


    // 2) compute the total size of all motifs
    for (const auto & motif:this->list_cms)
    {
        size_all_colomn_name+=motif.size()+1; // +1 for the separtor
    }
    size_all_colomn_name +=1; // for the endline of column name

    size_buffer+=size_all_colomn_name;

    // reserve the space memory , then add the informations
    buffer.reserve(size_buffer);

    // the column names:
    buffer+="familyId";
    for (const auto & motif:this->list_cms)
    {
        buffer+=csv_separator;
        buffer+=motif;
    }
    buffer+="\n";

    // the matrix values
    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->gst.isSequencesAreGroupedByFamilies())
        {
            auto pair_familyId_SeqId = this->get_FamilyId_And_SeqId(i); //(idx_family,idx_seq_in_family)
            buffer+=util::to_string(pair_familyId_SeqId.first); //idx_family
        }
        else
        {
            buffer+="0"; //idx_family is 0, we have one family
        }

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            buffer+=csv_separator;
            buffer+=util::to_string(matrix_nbOcrrs_cmsSeqsIds[i][j]); // nb occurence
        }
        buffer+="\n";
    }

    //-------------------------------------------------------------------

    // Open the file
    file_csv.open(file_output+".csv", ios::out | ios::binary);

    cout<<" file_csv.open (\"matrix_motifs_seqIds.csv\"); "<<endl;

    // Write the buffer to the file
    file_csv.write(buffer.c_str(), buffer.size());

    file_csv.close();
    cout<<"file_csv.close();"<<endl;
}


// chatGPT: save dirctly to the file stream
void CommonMotifs::saveMatrixCMS_ToCsv_File_dircrly(string file_output) const
{
    unsigned int nb_seqs= this->nb_all_sequences;
    unsigned int nb_motifs = this->list_cms.size();

    cout<<" Total nb_seqs = "<<nb_seqs<<endl;
    cout<<" Total nb_motifs = "<<nb_motifs<<endl;

    ofstream file_csv;
    const string csv_separator =",";
    
    // Open the file
    file_csv.open(file_output, ios::out);


    // the column names: // Write the header row
    file_csv<<"familyId";
    for (const auto & motif:this->list_cms)
    {
        file_csv<< csv_separator << motif;
    }
    file_csv << endl;

    // the matrix values : // Write the matrix values one row at a time
    for (unsigned int i=0; i<nb_seqs; ++i)
    {
        if(this->gst.isSequencesAreGroupedByFamilies())
        {
            auto pair_familyId_SeqId = this->get_FamilyId_And_SeqId(i); //(idx_family,idx_seq_in_family)
            file_csv << pair_familyId_SeqId.first; //idx_family, // chatGPT std::ofstream has an operator << that accepts integers as well as other data types and it writes them to the file as a string, it internally converts the integer to string before it writes to the file. So it's not an error, it's a built-in feature of std::ofstream.
        }
        else
        {
            file_csv << "0"; //idx_family is 0, we have one family
        }

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            file_csv << csv_separator << matrix_nbOcrrs_cmsSeqsIds[i][j]; // nb occurence
        }
        file_csv << endl;
    }

    //-------------------------------------------------------------------

    
    // Close the file
    file_csv.close();
    cout<<"file_csv.close();"<<endl;
}


// chatGPT: write by row, each row go to a buffer, then use stream write
// take the code generated by chatGPT directly
void CommonMotifs::saveMatrixCMS_ToCsv_File_rowByrow(string file_output) const
{
    unsigned int nb_seqs = this->nb_all_sequences;
    unsigned int nb_motifs = this->list_cms.size();

    cout<<" Total nb_seqs = "<<nb_seqs<<endl;
    cout<<" Total nb_motifs = "<<nb_motifs<<endl;

    std::ofstream myFile;
    myFile.open(file_output+".csv");
    // Write the header row
    std::string header;
    header.reserve(nb_motifs*20 + 10);
    header += "familyId";
    for (const auto & motif:this->list_cms)
    {
        header += ",";
        header += motif;
    }
    header += '\n';
    myFile.write(header.c_str(), header.size());

    std::string buffer;
    buffer.reserve(nb_motifs*5 + 10); // reserve space for the row
    // Write the matrix values
    for (unsigned int i = 0; i < nb_seqs; ++i)
    {
        buffer.clear();
        if(this->gst.isSequencesAreGroupedByFamilies())
        {
            auto pair_familyId_SeqId = this->get_FamilyId_And_SeqId(i); //(idx_family,idx_seq_in_family)
            buffer+=std::to_string(pair_familyId_SeqId.first); //idx_family
        }
        else
        {
            buffer+="0"; //idx_family is 0, we have one family
        }

        for (unsigned int j = 0; j < nb_motifs; ++j)
        {
            buffer+=",";
            buffer+=std::to_string(matrix_nbOcrrs_cmsSeqsIds[i][j]); // nb occurence
        }
        buffer+='\n';
        myFile.write(buffer.c_str(), buffer.size());
    }
    myFile.close();
    cout<<"file_csv.close();"<<endl;
}

// generated by chatGPT
// save matrix by chunks
// todo: check this, I am not sure that is correct.
void CommonMotifs::saveMatrixCMS_ToCsv_File_Chunked(string file_output, unsigned int chunk_size) const {
    unsigned int nb_seqs = this->nb_all_sequences;
    unsigned int nb_motifs = this->list_cms.size();

    cout << " Total nb_seqs = " << nb_seqs << endl;
    cout << " Total nb_motifs = " << nb_motifs << endl;

    std::ofstream myFile;
    myFile.open(file_output + ".csv");
    // Write the header row
    std::string header;
    header.reserve(nb_motifs * 20 + 10);
    header += "familyId";
    for (const auto &motif : this->list_cms) {
        header += ",";
        header += motif;
    }
    header += '\n';
    myFile.write(header.c_str(), header.size());

    std::string buffer;
    buffer.reserve(nb_motifs * 5 + 10); // reserve space for the row

    // Write the matrix values
    for (unsigned int i = 0; i < nb_seqs; i += chunk_size) {
        buffer.clear();
        for (unsigned int j = 0; j < chunk_size && (i + j) < nb_seqs; ++j) {
            if (this->gst.isSequencesAreGroupedByFamilies()) {
                auto pair_familyId_SeqId = this->get_FamilyId_And_SeqId(i + j); //(idx_family,idx_seq_in_family)
                buffer += std::to_string(pair_familyId_SeqId.first); //idx_family
            } else {
                buffer += "0"; //idx_family is 0, we have one family
            }

            for (unsigned int k = 0; k < nb_motifs; ++k) {
                buffer += ",";
                buffer += std::to_string(matrix_nbOcrrs_cmsSeqsIds[i + j][k]); // nb occurence
            }
            buffer += '\n';
        }
        myFile.write(buffer.c_str(), buffer.size());
    }
    myFile.close();
    cout << "file_csv.close();" << endl;
}


const unsigned int CommonMotifs::get_nb_totals_cms() const
{
    return this->list_cms.size();
}

// the total number of sequences 
// is the same thing as 
// the number of row in the matrix_nbOcrrs_cmsSeqsIds
unsigned int CommonMotifs::get_nb_all_sequences() const
{
    return this->nb_all_sequences;

    // number of row in the matrix_nbOcrrs_cmsSeqsIds
    // return this->matrix_nbOcrrs_cmsSeqsIds.size();
}

pair<unsigned int, unsigned int> CommonMotifs::get_FamilyId_And_SeqId(unsigned int global_seq_id) const
{   
    return SequenceIdManager::map_globalSeqId_To_FamilyAndLocalIds_Incremental_IndexBased(this->list_sum_nb_seqs, global_seq_id);
}

// ----------------------------------------------------------------------------
// more save to csv functions
// ----------------------------------------------------------------------------

void CommonMotifs::saveMatrixCMS_ToCsv_File_Optimized(const std::string& file_output) const {
    const unsigned int nb_seqs = this->nb_all_sequences;
    const unsigned int nb_motifs = this->list_cms.size();
    const std::string csv_separator = ",";

    std::ofstream file_csv(file_output, std::ios::out | std::ios::binary);
    if (!file_csv) {
        throw std::runtime_error("Unable to open file for writing: " + file_output);
    }

    // Write header
    file_csv << "familyId";
    for (const auto& motif : this->list_cms) {
        file_csv << csv_separator << motif;
    }
    file_csv << '\n';

    // Prepare buffer for rows
    const size_t row_buffer_size = nb_motifs * 10 + 100; // Estimate size
    std::string row_buffer;
    row_buffer.reserve(row_buffer_size);

    // Write matrix rows
    for (unsigned int i = 0; i < nb_seqs; ++i) {
        row_buffer.clear();

        // Add family ID
        ///* // no need to check.
        if (this->gst.isSequencesAreGroupedByFamilies()) {
            auto pair_familyId_SeqId = this->get_FamilyId_And_SeqId(i);
            row_buffer += std::to_string(pair_familyId_SeqId.first);
        } else {
            row_buffer += '0'; // idx_family is 0, we have one family
        }
        //*/
       //auto pair_familyId_SeqId = this->get_FamilyId_And_SeqId(i);
       //row_buffer += std::to_string(pair_familyId_SeqId.first);

        // Add motif occurrences
        for (unsigned int j = 0; j < nb_motifs; ++j) {
            row_buffer += csv_separator;
            row_buffer += std::to_string(matrix_nbOcrrs_cmsSeqsIds[i][j]);
        }
        row_buffer += '\n';

        file_csv.write(row_buffer.c_str(), row_buffer.size());
    }

    file_csv.close();
}

