//
// Created by ibra on 11/20/2019.
//

#include <algorithm>
#include <functional>
#include "SuffixTree_QuadraticTime.h"
#include "SequenceIdManager.h"

SuffixTree_QuadraticTime::SuffixTree_QuadraticTime()
{
    help_initializer();

    this->max_length_motif=0; // this mean we don't use max, we use for each sequences its length.
    this->min_length_motif=1; // the default. only one char

    cout<<"default constructor"<<endl;
}

SuffixTree_QuadraticTime::SuffixTree_QuadraticTime(uint32_t maxLengthMotif, uint32_t minLengthMotif) : max_length_motif(
        maxLengthMotif)
{
    if (minLengthMotif==0) this->min_length_motif=1; // the min default is one char.
    else
        this->min_length_motif = minLengthMotif;

    if(maxLengthMotif!=0 && maxLengthMotif < minLengthMotif ) this->max_length_motif = minLengthMotif; // its a choice.

    help_initializer();
}

void SuffixTree_QuadraticTime::help_initializer() //help initializer
{
    this->root_suffix_tree=new Node('\0');
    this->are_sequences_grouped_in_families = false;
    this->nb_all_sequences=0;
}

SuffixTree_QuadraticTime::~SuffixTree_QuadraticTime() {

    //cout<<"delete is called on SuffixTree_QuadraticTime"<<endl;
    if(this->root_suffix_tree) // here we must to check because in code of Node we acces to childer whithout cheking the nullptr
    {
        delete(this->root_suffix_tree);
    }

    list_seqs.clear(); // 1) all sequences independent from initial clusturing
    list_families_seqs.clear();     // 2) sequnces are grouped in initial families (cluster)
    list_sum_nb_seqs.clear(); // for indexing seqs by family
}

/**
 *
 * @param str  the string to add to the tree (simple suffixe tree or genrilized )
 * @param idx_seq
 */
void SuffixTree_QuadraticTime::AddString(string str, unsigned int idx_seq)
{
    //cout<<"AddString: "<<str<<endl;

    int num_char=-1; // [0..51]<==>[A-Za-z]

    Node* my_node= nullptr;


    // ajouter le mot à l'arbre caracter par caracter.
    for(unsigned int i=0;i<str.size();i++)
    {
        //cout<<"insert : ";

        my_node=this->root_suffix_tree;

        for (unsigned int j = i; j <str.size() ; ++j)
        {
            //cout<<str[j];

            num_char = get_num_char(str[j]);
            if(num_char==-1) continue;

            // si le char n'existe pas parmet les fills dc add new node
            if( my_node->listChildren[ num_char ]== nullptr)
            {
                my_node->listChildren[ num_char ] = new Node(str[j]);
                my_node=my_node->listChildren[ num_char ];
            }
            else // donc on a un fils qui a ce char, donc on recupaire le node fils de cette branche et continue le traitement
            {
                my_node=my_node->listChildren[ num_char ];
            }
        }

        my_node->AddNewPosition(idx_seq, i);

        //cout<<endl;
    }

    ///my_node->nb_word = pos_word; // sauvgarder la position de mot (sa position dans tab_dic)
}

void SuffixTree_QuadraticTime::GenerateSuffixTree(string &str)
{
    //cout<<"GenerateSuffixTree: "<<str<<endl;
    //this->seq=str; // in this code we don't use this.
    //this->AddString(str, 0);
    this->addAllSubMotifMinMax(str, 0);
}

/**
 *
 * we use only 4 chars (ACTG)
 * for the others combinaison (UPAC) codes, on va pas les represnter
 * on va par exemple si on trouve un "N" donc tous les 4 caracters, on va generer pour les 4 fils en meme temps
 * sachant que ces caractere il y a pas beacoup, donc comme cela on gagne de l'espace
 * c'est mieu que de faire pour chaque node un tab de 15 fils ou plus.
 *
 * @param my_char
 * @return
 */
// TODO: use a map instead of switch case
int SuffixTree_QuadraticTime::get_num_char(char my_char) const
{
    switch  (my_char)
    {
        case 'A' : return  0;
        case 'a' : return  0;
        case 'T' : return  1;
        case 't' : return  1;
        case 'U' : return  1;
        case 'u' : return  1;
        case 'G' : return  2;
        case 'g' : return  2;
        case 'C' : return  3;
        case 'c' : return  3;
            /*
            case 'S' : return  4;
            case 'W' : return  5;
            case 'R' : return  6;
            case 'Y' : return  7;
            case 'K' : return  8;
            case 'M' : return  9;
            case 'B' : return 10;
            case 'V' : return 11;
            case 'H' : return 12;
            case 'D' : return 13;
            case 'N' : return 14;
            case '-' : return 15;
             */
        default  : return -1;
    }
}

const Node *SuffixTree_QuadraticTime::getRootSuffixTree() const {
    return root_suffix_tree;
}

bool SuffixTree_QuadraticTime::CheckSuffixTree(string str) const
{
        cout << "CheckSuffixTree: " << str << endl;

        int num_char; // [0..51]<==>[A-Za-z]

        Node *my_node = nullptr;


        // ajouter le mot à l'arbre caracter par caracter.
        for (unsigned int i = 0; i < str.size(); i++)
        {
            cout << "check : ";

            my_node = this->root_suffix_tree;

            for (unsigned int j = i; j < str.size(); ++j)
            {
                cout << str[j];

                num_char = get_num_char(str[j]);

                // si le char n'existe pas parmet les fills dc add new node
                if (my_node->listChildren[num_char] == nullptr)
                {
                    return false;
                } else // donc on a un fils qui a ce char, donc on recupaire le node fils de cette branche et continue le traitement
                {
                    my_node = my_node->listChildren[num_char];
                }
            }

            cout << endl;
        }

        ///my_node->nb_word = pos_word; // sauvgarder la position de mot (sa position dans tab_dic)

        cout<<"the string "<<str<<" exist "<<endl;
        return true;
}

/*
    Generate a suffix tree for a list of strings.
    one list of strings, they are not grouped in families.
*/
void SuffixTree_QuadraticTime::GenerateGeneralizedSuffixTree(vector<string> &list_strings)
{
    this->list_seqs=list_strings;
    this->nb_all_sequences= list_strings.size();
    this->are_sequences_grouped_in_families = false; // since we have only one list of strings

    for (unsigned int i = 0; i <list_strings.size() ; ++i)
    {
        //this->AddString(list_strings.at(i), i);
        this->addAllSubMotifMinMax(list_strings.at(i), i);
    }

}

/**
 * Generate a generalized suffix tree for a list of strings, grouped in families.
 *  */
void SuffixTree_QuadraticTime::GenerateGeneralizedSuffixTree(vector<vector<string>> &list_grouped_strings)
{
    this->list_families_seqs=list_grouped_strings;
    this->are_sequences_grouped_in_families = true;

    unsigned int global_sequence_id=0; // global sequence id = index of the seq in all sequences of all families. just incrementing by one for each sequence

    for (const auto & list_single_family_seqs : list_families_seqs)
    {
        for (const auto & sequence : list_single_family_seqs)
        {
            //this->AddString(sequence, sequence_id);
            this->addAllSubMotifMinMax(sequence, global_sequence_id);
            global_sequence_id++;
        }

        nb_all_sequences+=list_single_family_seqs.size();
    }
    
    this->build_family_seq_global_index_map();// ibra, we can integrate this function here,
                                // because they have the same loop, sauf si on vaut laisser la visibilite.


    cout << "nb_all_sequences after GenerateGeneralizedSuffixTree: " << nb_all_sequences << endl;
    cout << "list_families_seqs.size(): " << list_families_seqs.size() << endl;
    cout << "list_sum_nb_seqs.size(): " << list_sum_nb_seqs.size() << endl;
    cout << "list_sum_nb_seqs: ";
    for (const auto &sum : list_sum_nb_seqs) {
        cout << sum << " ";
    }
    cout << endl;
}

void SuffixTree_QuadraticTime::build_family_seq_global_index_map()
{
    vector<unsigned int> list_nb_elements_per_family;
    for (const auto& family : list_families_seqs) {
        list_nb_elements_per_family.push_back(family.size());
    }
    this->list_sum_nb_seqs = SequenceIdManager::generate_list_sum_nb_elements(list_nb_elements_per_family);
}

bool SuffixTree_QuadraticTime::CheckGST(vector<string> &listStrings) const
{
    for(const auto &str:listStrings)
    {
        if(!this->CheckSuffixTree(str)) return false;
    }

    cout<<"all the strings exist in the GST yopyy "<<endl;

    return true;
}

const vector<string> &SuffixTree_QuadraticTime::getListSeqs() const {
    return list_seqs;
}

const vector<vector<string>> &SuffixTree_QuadraticTime::getListFamiliesSeqs() const {
    return list_families_seqs;
}



/**
 *
 * @param id_seq the global id of the sequence
 * @return a reference to the location of sequence (where it is saved)
 */
const string &SuffixTree_QuadraticTime::getSeq(unsigned int global_seq_id) const
{
    // if we use the only sequences without families
    if(!this->are_sequences_grouped_in_families)
        return this->list_seqs.at(global_seq_id);

    // else, sequences are grouped in families
    pair<unsigned int,unsigned int> ij = this->get_FamilyId_And_SeqId(global_seq_id);

    return this->list_families_seqs.at(ij.first).at(ij.second);
}

bool SuffixTree_QuadraticTime::isSequencesAreGroupedByFamilies() const {
    return are_sequences_grouped_in_families;
}



//pair<unsigned int, unsigned int> SuffixTree_QuadraticTime::get_FamilyId_And_SeqId(unsigned int seq_id) const
pair<unsigned int, unsigned int> SuffixTree_QuadraticTime::get_FamilyId_And_SeqId(unsigned int global_seq_id) const
{
    // here list_sum_nb_seqs is empty, I don't know why ??? 
    std::cout << "in SuffixTree_QuadraticTime::get_FamilyId_And_SeqId" << std::endl;
    std::cout << "global_seq_id: " << global_seq_id << std::endl;
    std::cout << "list_sum_nb_seqs size: " << this->list_sum_nb_seqs.size() << std::endl;
    std::cout << "list_sum_nb_seqs: ";
    for (const auto &elem : this->list_sum_nb_seqs) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
    std::cout << " just befor call map_globalSeqId_To_FamilyAndLocalIds_Incremental_IndexBased" << std::endl;
    
    return SequenceIdManager::map_globalSeqId_To_FamilyAndLocalIds_Incremental_IndexBased(this->list_sum_nb_seqs, global_seq_id);
}


unsigned int SuffixTree_QuadraticTime::getNbAllSequences() const{
    return nb_all_sequences;
}

void SuffixTree_QuadraticTime::addAllSubMotifMinMax(const string &str, uint32_t idx_str)
{
    if (str.size() < this->min_length_motif) return;

//    cout<<"addAllSubMotifMinMax : "<<endl;
//    cout<<"str : "<<str<<endl;
//    cout<<"idx_str: "<<idx_str<<endl;
//    cout<<"min : "<<min_length_motif<<endl;
//    cout<<"max : "<<max_length_motif<<endl;

    // maxLengthMotif==0, mean we don't need max, so call diffirent method, to avoid to check each time when extracting motif
    if (this->max_length_motif == 0) return this->addAllSubMotifMin(str, idx_str);

    // it is maybe we find some strings that have size() < this->max_length_motif
    // for that we use a local_str_max_length_motif that have the size of those str, and we don't alter the global this->max_length_motif
    uint32_t local_str_max_length_motif = this->max_length_motif;
    if(str.size() < local_str_max_length_motif) local_str_max_length_motif = str.size();

    uint32_t last = str.size() - this->min_length_motif + 1;

    uint32_t last_max = str.size() - local_str_max_length_motif + 1;

    for (uint32_t i = 0; i < last_max; ++i)
    {
        addSubStr(str, idx_str, i, i + local_str_max_length_motif);
    }

    for (uint32_t i = last_max; i < last; ++i)
    {
        addSubStr(str, idx_str, i, str.size());
    }
}

void SuffixTree_QuadraticTime::addAllSubMotifMin(const string &str, uint32_t idx_str)
{
    uint32_t str_size = str.size();
    uint32_t last_pos_min_motif = str.size() - this->min_length_motif + 1;

    for (uint32_t i = 0; i < last_pos_min_motif; ++i)
    {
        addSubStr(str, idx_str, i, str_size);
    }
}

void SuffixTree_QuadraticTime::addSubStr(const string &str, uint32_t idx_str, uint32_t from, uint32_t to)
{
    // cout<<" addSubStr : "<<endl;
    // cout<<"str : "<<str.substr(from,to-from)<<endl;
    // cout<<"idx_str: "<<idx_str<<endl;
    // cout<<"char by char insert : ";

    int num_char; // [0..51]<==>[A-Za-z]

    Node* my_node = this->root_suffix_tree;

    // add the subStr[from--to[ to generilized suffixe tree.
    for(uint32_t i = from ; i < to ; ++i)
    {
        // cout<<str[i];

        num_char = get_num_char(str[i]);
        if(num_char==-1) continue; // to not consider the other wrong chars, we just ignore it and jump to to the next char.

        // si le char n'existe pas parmet les fills dc add new node
        if(my_node->listChildren[ num_char ] == nullptr)
        {
            my_node->listChildren[ num_char ] = new Node(str[i]);
            my_node=my_node->listChildren[ num_char ];
        }
        else // donc on a un fils qui a ce char, donc on recupaire le node fils de cette branche et continue le traitement
        {
            my_node=my_node->listChildren[ num_char ];
        }

    }
    //cout<<endl;

    my_node->AddNewPosition(idx_str, from);
}

uint32_t SuffixTree_QuadraticTime::getMinLengthMotif() const {
    return min_length_motif;
}

uint32_t SuffixTree_QuadraticTime::getMaxLengthMotif() const {
    return this->max_length_motif;
}

const vector<unsigned int> SuffixTree_QuadraticTime::getListSumNbSeqs() const {
    return list_sum_nb_seqs;
}
