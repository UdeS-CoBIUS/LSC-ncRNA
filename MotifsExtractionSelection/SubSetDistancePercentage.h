//
// Created by ibra on 1/21/2020.
//

#ifndef TEST_SUBSETDISTANCEPERCENTAGE_H
#define TEST_SUBSETDISTANCEPERCENTAGE_H

#include <unordered_map>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//
// problem definition
// let be a list of n no negative integers, X={x1, x2,...,xn}.
// let alpha be a threshold to check if abs(x_i - x_j)<alpha (with i=!j and i,j \in [1..n])
// let a subset Y={y1,y2,...,ym} with Y in X (m<=n, element of Y are in X)
// Task 1:
// --------
// find a subset Y (all subset Y) where all its elements respect the threshold alpha
// abs(y_i-y_j)<alpha (y_i,yj \in [1..m] and i!=j).
//     All elements in X that respect the alpha threshold must be gathered in Y.
//
// let beta:  a percentage (a threshold), to check if the count(Y)=m (number of element) is acceptable comparing to count(X)=n
// (m*100)/n >= beta
//
// task 2:
// -------
// find a set Y that respect the threshold alpha and beta.
//
// task 3:
// in the case where there are multiple subset Y, find the best subset
// the best subset!! ==
//          1) subset that have the bigest 'm' (count(Y))
//          2) the subset that have the gretest element x (x \in X)
//             or the best centroid "c" where c is the moyenne in Y, or [c-alpha<--c-->c+alpha]
//              because each subset Y, contain only element that are in the range [c-alpha<--c-->c+alpha] with abs(maxvalueInRange - minValueInRange)<alpha
//
// ------------------------------------------------
// example 1:
// X={5,6,5,4,8,9,10,9,1,9} ; count(X)=10=n
// alpha=1
// beta=40%
//
// task1: find all Y that respect alpha.
//  >> each element alone in X can constitue a subset Y,
//     because one element can't be checked against alpha.
//     It is a choice, we can and take it as subset or not.
//
//    the subsets can't be taken excepte Y_0_9:{1},
//    because they can be merged with other subsets and still respect the threshold aplha.
//    Y_0_1:{5},Y_0_2:{6},Y_0_3:{5},Y_0_4:{4},Y_0_5:{8},Y_0_6:{9},
//    Y_0_7:{10},Y_0_8:{9},Y_0_9:{1},Y_0_10:{9}
//
//  Results:
// Y0={1}
// Y1={5,5,6}
// Y2={5,5,4}     Y1 and Y2 can't be one subset because 6-4=2 > alpha =1
// Y3={8,9,9,9}
// Y4={9,9,9,10}
//
// task 2: choose a subset that respect beta.  (>=beta)
// Y3={8,9,9,9}
// Y4={9,9,9,10}
//
// task 3: choose the best Y
// i think Y4={9,9,9,10} that have moyenne > Y3
//   but the centroidn in both cases i think is 9.
// we can see that clearly if alpha = 2.
//
// ------------------------------------------------
// example 2:  same as 1, execpt we change alpha
// X={5,6,5,4,8,9,10,9,1,9} ; count(X)=10=n
// alpha=2
// beta=40%
//
// task 1: all Y that respect alpha
// Y0={1}
// Y1={5,5,6,4}
// Y2={6,8}
// Y3={8,9,9,9,10}
//
// task 2: choose a subset that respect beta.  (>=beta)
// Y1={5,5,6,4}
// Y3={8,9,9,9,10}
//
// task 3: choose the best Y
// Y3={8,9,9,9,10}
//
//
// ibra solution in O(alpha*n) hihihi:
// --------------
//
// how to construct one Y:
// 1) To form the subset Y, we begin by the first element x_i.
// 2) this imply the possible elements that can be added to Y are in the range
//    [alpha-x_i,x_i+alpha]
// 3) we can have in Y tow element x_i+alpha and x_i-alpha, because we violate the condition
//    so, the maxValieInRange-minValueInRange <=alpha
//   before addin an new element to Y, we have to ensure that it do not violate the condition
//   for that, wa can for example check the x_j that we want to insert with all elemenst in Y.
//  Or better:
//  before adding a new ellement x_j, compar it only with the minValue and maxValue
//  how: .....
//      initialize min and max by the x_i (wich is the key)
//      check if (abs(x_j-max)<=alpha && (abs(x_j-min)<=alpha) add it to Y.
//      update min, max: if we can add x_j to Y,
//                       we check if x_j>max, max = x_j
//                       else
//                       we check if x_j<min, min = x_j
//         je ne sais pas si je peut faire cela mieux, reduire le nombre des verification...
// 4) loop thought all element in X, to add them in Y.
//
// how to construct all the subset Y:
// for all element x_i in X:
// do: 1,2,3,4.
// this mean O(n^2)
//
/// now solution in O(alpha*n)
//
// using: map (hash Table (key, value))
// key : is the integer x_i in X.
// value: {a set Y{}}
// use a tmporary var: minValue, maxValue
//
// first setp:
// -----------
// initialize the map with enly the key. (the value remain empty)
// this take O(n), loop thought all the X elements.
// the repeated element in X will have the same entry in the hashTble.
// this is good, because abs(x_i - x_i) = 0 <= alpha , because alpha is alwys >=0.
//
//
//
// second step:
// ------------
// for each element x_i in X (loop thought all the X elements)
// in the range [x_i-alpha,x_i+alpha], we use y_j \in the range.
//                      we can restrect the range according to nown information about X
//                      example , if X contain only postivie number,
//                      so, if x_i-alpha<0, begin from 0.
// A) if y_j exist in the hashmap:
// B) do 3 in prevoius algorithm.
//
// -----------------------------------------------------------
// pire case example :
// X={0,1,2,3,...,n}
// we will have n-alpha subset Y.
// each time we have window thta slide from left to right that contain alpha+1 element.
//
// best case example:
// a set that contain the same element repeated n times.
// this mean we will have only 1 entry in the hashmap
// so as result, the set Y=X. (independament de changement de alpha et beta)
// we we can check if the siz of hashmap = 1, return X.
// this mean we just stop at end of the initialization.
//
// otherwise: pire case, best case, we have the sampe number of operation. (i think)
//

class SubSetDistancePercentage{

    struct SubSet{
        int min=0;
        int max=0;
        vector<int> list;

        SubSet() = default;

        SubSet(int min, int max) : min(min), max(max) {}
    };

    template < typename T > static std::string to_string(const vector<T> &vec)
    {
        //if(&vec== nullptr || vec.empty()) return "{}";
        if(vec.empty()) return "{}";

        std::ostringstream stm ;
        stm<<"{";

        for(uint32_t idx=0; idx<vec.size()-1;idx++)
        {
            stm<<vec.at(idx);
            stm<<", ";
        }

        stm<<vec.at(vec.size()-1);

        stm<<"}";

        return stm.str();
    }

public:


// unsigned int x,y;
// x-y : bad idea
// abs(x-y), also bad ida.
// https://stackoverflow.com/questions/27833289/does-absunsigned-long-make-any-sense

/**
 * tell f we can add val or not, and update minMax value.
 * @param min the minValue in that exist in the set Y
 * @param max the maxValue in that exist in the set Y
 * @param alpha a threshold to accept new value in Y
 * @param keyValue the first value add to Y. and it is a key in hashmap, or the centre of all value add after [alpha-keyVal,alpha+keyValue]
 * @param val  the new value that we want to add
 * @return can we add val or not.
 */
    //static inline bool canAddValTolist_andUpdateMinMax(uint32_t &min, uint32_t &max, uint32_t alpha, uint32_t keyValue, uint32_t val)
    static inline bool canAddValTolist_andUpdateMinMax(int &min, int &max, int alpha, int keyValue, int val)
    {

        if(val < keyValue) // ...min...val...KeyValue...max...> // if max==KeyValue, its remain correct, because KeyValue-val, true if val is between min and Keyvalue, true if val==min, true if val<min and KeyValue-val<=alpha
        {
            if(max-val<=alpha)
            {
                if(val<min) min = val;
                return true;
            }
        }
        else // ...min...KeyValue...val...max...>
        {
            if (val - min <= alpha)
            {
                if (val > max) max = val;
                return true;
            }
        }

        return false;
    }


    /**
     *
     * @param list_integers
     * @param alpha
     *        beta , it is not used here, it is used by the function that call this method
     * @return a pair of best_key (centroid), list best sub set acording to alpha and beta
     */
    static pair<int,vector<int>> getBestSubSet(const vector<int> &list_integers ,int alpha)
    {

        int nb_all_elements = list_integers.size();

        unordered_map<int, SubSet > map_keyValue_listRespectAlpha; // keyvalue ; {(min,max),liste_element}

        for (int key_value : list_integers)
        {
            if(map_keyValue_listRespectAlpha.find(key_value) == map_keyValue_listRespectAlpha.end())
            {
                map_keyValue_listRespectAlpha.insert(make_pair(key_value,SubSet(key_value,key_value) )); // we insert with 0, because, just after we will add+1 in the range of [+-alpha]
            }
        }


        for (int val : list_integers)
        {
            int begin=0;
            if(val - alpha > begin) begin= val - alpha;
            int end = val + alpha;


            for (int key_value = begin; key_value <= end; ++key_value)
            {
                // check if val exist as element in our list
                if(map_keyValue_listRespectAlpha.find(key_value) != map_keyValue_listRespectAlpha.end())
                {
                    if(canAddValTolist_andUpdateMinMax(map_keyValue_listRespectAlpha[key_value].min,map_keyValue_listRespectAlpha[key_value].max,alpha,key_value,val))
                    {
                        map_keyValue_listRespectAlpha[key_value].list.push_back(val);
                    }
                }
            }
        }


        // get the elemet of max score
        int best_keyValue=0; // we initilize with worst keyValue, in our case here is 0
        int best_keyValue_subSet_size=0;
        vector<int > best_list;
        for (auto &oneMapElment : map_keyValue_listRespectAlpha)
        {
            if (oneMapElment.second.list.size() >best_keyValue_subSet_size)
            {
                best_keyValue = oneMapElment.first;
                best_list = oneMapElment.second.list;
                best_keyValue_subSet_size = oneMapElment.second.list.size();
            }
            else if(oneMapElment.second.list.size() == best_keyValue_subSet_size && oneMapElment.first>best_keyValue )
            {
                best_keyValue = oneMapElment.first;
                best_list = oneMapElment.second.list;
            }
        }

        return make_pair(best_keyValue,best_list);
    }

    // use another paramater : the lower bound allowed
    static pair<int,vector<int>> getBestSubSet_lowerBound(const vector<int> &list_integers ,int alpha , unsigned int lower_bound_allowed)
    {

        int nb_all_elements = list_integers.size();

        unordered_map<int, SubSet > map_keyValue_listRespectAlpha; // keyvalue ; {(min,max),liste_element}


        for (int key_value : list_integers)
        {
            // skip the key_value that are < lower_bound_allowed
            if(key_value >= lower_bound_allowed &&
               map_keyValue_listRespectAlpha.find(key_value) == map_keyValue_listRespectAlpha.end())
            {
                map_keyValue_listRespectAlpha.insert(make_pair(key_value,SubSet(key_value,key_value) )); // we insert with 0, because, just after we will add+1 in the range of [+-alpha]
            }
        }

        for (int val : list_integers)
        {
            if(val < lower_bound_allowed) continue;  // skip the val that are < lower_bound_allowed

            int begin=0;
            if(val - alpha > begin) begin= val - alpha;
            int end = val + alpha;

            for (int key_value = begin; key_value <= end; ++key_value)
            {
                // check if val exist as element in our list
                if(map_keyValue_listRespectAlpha.find(key_value) != map_keyValue_listRespectAlpha.end())
                {
                    if(canAddValTolist_andUpdateMinMax(map_keyValue_listRespectAlpha[key_value].min,map_keyValue_listRespectAlpha[key_value].max,alpha,key_value,val))
                    {
                        map_keyValue_listRespectAlpha[key_value].list.push_back(val);
                    }
                }
            }
        }


        // Best: based on the length of the set, bigest length is the best
        int best_keyValue=0; // we initilize with worst keyValue, in our case here is 0
        int best_keyValue_subSet_size=0;
        vector<int > best_list;
        for (auto &oneMapElment : map_keyValue_listRespectAlpha)
        {
            if (oneMapElment.second.list.size() >best_keyValue_subSet_size)
            {
                best_keyValue = oneMapElment.first;
                best_list = oneMapElment.second.list;
                best_keyValue_subSet_size = oneMapElment.second.list.size();
            }
            else if(oneMapElment.second.list.size() == best_keyValue_subSet_size && oneMapElment.first>best_keyValue )
            {
                best_keyValue = oneMapElment.first;
                best_list = oneMapElment.second.list;
            }
        }

        return make_pair(best_keyValue,best_list);
    }

};

#endif //TEST_SUBSETDISTANCEPERCENTAGE_H
