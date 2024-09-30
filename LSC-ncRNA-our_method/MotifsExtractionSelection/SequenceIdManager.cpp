#include "SequenceIdManager.h"


// implementation of the methods

/**
 * given a seq id, access the sequence in vector<vector<string>>.
 * where all seqs have id from 0 to n (n is the total number of all seqs in all families)
 * how we can acheave this:
 * index the sequences id by families (vector),
 * by putting the last seq id of each vector in bucket, from the first to the last vector
 * after we use getFamilyId to retrive the index of the vector acoording to seq id
 * we use decotomy search (the function lower_bound).
 *
 * we need -1 at first, since it will be used as index for seq by family, the size give the size number of seqs, but we begin from 0, so last seq is size -1 (for 1st family), then for others, we just add size. after given family id, and local seq id in family, we can retrive global index using list_sum
 *
 * @param list_nb_elements_per_family: list of the number of elements in each family
 * @return list of the sum of the number of elements in each family
 */
vector<unsigned int> SequenceIdManager::generate_list_sum_nb_elements(const vector<unsigned int> &list_nb_elements_per_family)
{
    vector<unsigned int> list_sum_nb_elements;
    list_sum_nb_elements.reserve(list_nb_elements_per_family.size()); // Reserve space.

    unsigned int sum = list_nb_elements_per_family[0] - 1; // first  family size. -1 because we begin from 0 for indexing.
    list_sum_nb_elements.push_back(sum);

    for (size_t i = 1; i < list_nb_elements_per_family.size(); ++i)
    {
        sum += list_nb_elements_per_family[i];
        list_sum_nb_elements.push_back(sum);
    }

    return list_sum_nb_elements;
}


/**
 *
 * @param global_seq_id global sequence id
 * @return the position of the sequence in wich families (the index of vector) in vector<vector<string>>
 * and its index in its own local family, which mean in vector<string>
 */
pair<unsigned int, unsigned int>
SequenceIdManager::map_globalSeqId_To_FamilyAndLocalIds_Incremental_IndexBased(const vector<unsigned int> &list_sum_nb_elements,
                                                            unsigned int global_seq_id)
{
    //std::cout << "we are in get_FamilyId_And_SeqId_IndexBased" << std::endl;
    //std::cout << "seq_id in: " << seq_id << std::endl;

    unsigned int idx_family=0; // idx vector
    unsigned int idx_seq_in_family=0;

    // we dont't need to check the result of lower_bound, because all global_seq_id are (must be) valid, they exist in the range of list_sum_nb_seq
    idx_family = std::lower_bound(list_sum_nb_elements.begin(),list_sum_nb_elements.end(),global_seq_id)
                    - list_sum_nb_elements.begin();

    if(idx_family==0) idx_seq_in_family=global_seq_id;
    else{
        idx_seq_in_family = global_seq_id-list_sum_nb_elements.at(idx_family-1)-1;
    }

    // Print the result for debugging
    //std::cout << "idx_family: " << idx_family << std::endl;
    //std::cout << "idx_seq_in_family: " << idx_seq_in_family << std::endl;

    return make_pair(idx_family,idx_seq_in_family);
}

/**
 * return a id for sequence composed by idx_family and idx_seq_in_family
 * the seq id is unsigned int which is composed by tow part:
 * the first 16 bits stor the idx_family
 * the second 16 bits stor the idx_seq_in_family
 *
 * important: we use this method when:
 * 1) the number of sequence in each family <= 65,536 (2^16)
 * 2) the number of all family <= 65,536 (2^16)
 *
 * in practice, for now, these number is very suficient to holde the information
 * because, the number of sequence or family is << 65,536
 * for now Rfam 14.1 (January 2019, 3016 families)
 *
 * ofcourse, we supose the use of 32 bits for int. which the case for machine of now days.
 *
 * @param idx_family the id of the family, or the vector where the sequence is stored
 * @param idx_seq_in_family the id of the seq in the family (in the vector)
 * @return a unique id for seq
 */ 
unsigned int SequenceIdManager::compose_binary_global_seq_id(unsigned int idx_family, unsigned int idx_seq_in_family)
{
   // i will not cast the value for (uint32_t), because all variable are: unsigned int
    // and by definition uint_32 in stdint.h is (typedef unsigned		uint32_t;)
    // the shift by 16 must be unsigned hence 16u
    // https://stackoverflow.com/questions/1294649/cleanest-way-to-combine-two-shorts-to-an-int
    // https://stackoverflow.com/questions/50399090/use-of-a-signed-integer-operand-with-a-binary-bitwise-operator-when-using-un

    //uint32_t seq_id = (uint32_t)((uint32_t)((uint32_t)idx_family<<16) | (uint32_t)idx_seq_in_family);
    //unsigned int seq_id = (idx_family<<16u) | idx_seq_in_family;
    //return seq_id;

    if (idx_family > 0xFFFF || idx_seq_in_family > 0xFFFF) {
        throw std::out_of_range("Family or sequence index exceeds 16-bit limit");
    }

    return (idx_family<<16u) | idx_seq_in_family;
}



std::pair<unsigned int, unsigned int>
SequenceIdManager::decompose_binary_global_seq_id_to_family_and_localSeqIds(unsigned int binary_global_seq_id)
{
    // this code comment for more explanation
    //unsigned int idx_family= (seq_id & 0xFFFF0000u)>>16u;
    //unsigned int idx_seq_in_family=(seq_id & 0x0000FFFFu);
    //return make_pair(idx_family,idx_seq_in_family);

    return make_pair((binary_global_seq_id & 0xFFFF0000u)>>16u, (binary_global_seq_id & 0x0000FFFFu));
}