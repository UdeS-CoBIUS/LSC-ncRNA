//
// Created by ibra on 11/27/2019.
//

#ifndef SUFFIX_TREE_UTILE_H
#define SUFFIX_TREE_UTILE_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

namespace util
{
#define ALPHA 2
#define BETA 60                     // BETA is pecentage of finding cm in the family
#define IS_CHECK_COMMON_MOTIF_ONLY 0 // 1 true , 0 false
#define IS_USE_ALPHA 1
#define IS_USE_BETA 0
#define BETA_2 20                   // BETA_2 is pecentage of finding cm acrosse diffirent families

#define GAMMA 2 // the minimum number of sequences that have a given Common motif (cm),
               // a cm is at least exist in tow sequences
               // so deafult value is GAMMA 2
               // this used when constrcating the firt matrix of common motif
               // Note: GAMMA is related to BETA,
               // if we have 10 seqs, and BETA = 50 % ==> the cm must exist in at leat 5 seqs
               // so in this cas BETA include GAMMA
               // but,
               // 1) if BETA = 0, we must check GAMMA
               // 2) there is a lot of motif that exist only in 1 seq, so
               //    the test with GAMMA is faster than the test with BETA
               //    (one comparaison vs compute of percentage)
               // if (list_motif >= GAMMA  && is_accepted_acording_to_Beta...)
               // here if the first condition is not true, we don't compute the percentage Beta...

#define CSV_SEPARATOR ","

const char kPathSeparator =
#ifdef _WIN32
            '\\';
#else
    '/';
#endif

    // it dosn't work, it call ( const T& n ) instead.
    template < typename T > std::string to_string(const vector<T> *vec)
    {
        cout<<"call to string avec pointeur"<<endl;
        return util::to_string(*vec);
    }

    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }

    template < typename T > std::string to_string(const vector<T> &vec)
    {
        //if(&vec== nullptr || vec.empty()) return "{}";
        if(!(&vec) || vec.empty()) return "{}";

        std::ostringstream stm ;
        stm<<"{";

        for(int idx=0; idx<vec.size()-1;idx++)
        {
            stm<<vec.at(idx);
            stm<<", ";
        }

        stm<<vec.at(vec.size()-1);

        stm<<"}";

        return stm.str();
    }

    template < typename T > std::string to_string(const vector<T*> &vec)
    {
        if(!(&vec) || vec.empty()) return "{}";

        std::ostringstream stm ;
        stm<<"{";

        for(int idx=0; idx<vec.size()-1;idx++)
        {
            stm<<(*vec.at(idx));
            stm<<", ";
        }

        stm<<(*vec.at(vec.size()-1));

        stm<<"}";

        return stm.str();
    }

    template < typename T > std::string to_string(T *array, int size)
    {
        if(size==0) return "{}";

        std::ostringstream stm ;
        stm<<"{";

        for(int idx=0; idx<size-1;idx++)
        {
            stm<<array[idx];
            stm<<", ";
        }

        stm<<array[size-1];

        stm<<"}";

        return stm.str();
    }



    // not tested yet.
    template < typename T > void VectorConcat(vector<T> &vec_dest, vector<T> &vec_src)
    {
        vec_dest.insert(vec_dest.end(), vec_src.begin(), vec_src.end());
    }

    // not tested yet.
    template < typename T > void VectorConcat(vector<T> &vec_dest, vector<T> &vec_src_1,  vector<T> &vec_src_2)
    {
        VectorConcat(vec_dest,vec_src_1);
        VectorConcat(vec_dest,vec_src_2);
    }
}


#endif //SUFFIX_TREE_UTILE_H
