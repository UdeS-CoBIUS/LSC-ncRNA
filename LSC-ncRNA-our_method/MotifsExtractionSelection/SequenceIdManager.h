#include <utility> // Add this at the top of the file
#include <vector>

using namespace std;

class SequenceIdManager {
public:

    // methods 1:
    // based on the list of the sum of the number of elements in each family,
    // global id => simply incremented id from 0 to nb all sequences - 1, in the order of the families.
    // list_sum_nb_elements => list of the sum of the number of elements in each family: 
    // list_sum_nb_elements[i] = sum of the number of elements in the families before the i-th family

    static vector<unsigned int> generate_list_sum_nb_elements(const vector<unsigned int> &list_nb_elements_per_family);

    static pair<unsigned int, unsigned int> map_globalSeqId_To_FamilyAndLocalIds_Incremental_IndexBased(const vector<unsigned int> &list_sum_nb_elements,
                                                            unsigned int global_seq_id);

    // methods 2: binary composition of the family and the local id
    // combination of family and local id: in one unique global id = binary composition in one value right and left of the binary representation of the family and the local id
   
    static unsigned int compose_binary_global_seq_id(unsigned int idx_family, unsigned int idx_seq_in_family);
    // based on the global id, it returns the family and the local id in the family
    static std::pair<unsigned int, unsigned int> decompose_binary_global_seq_id_to_family_and_localSeqIds(unsigned int binary_global_seq_id) noexcept;
};
