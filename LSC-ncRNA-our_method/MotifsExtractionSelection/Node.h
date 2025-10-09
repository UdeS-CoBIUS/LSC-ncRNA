//
// Created by ibra on 11/20/2019.
//

#ifndef SUFFIX_TREE_NODE_H
#define SUFFIX_TREE_NODE_H

#include <ostream>
#include <sstream>
#include <unordered_map>
// #include <bits/unique_ptr.h> // it is not used in this code, also it is not part of standard library, (can use memory instead). There's no need to replace it with <memory> since the code is not using smart pointers.

#include "Utile.h"

#define NB_CHAR 4
//#define END_CHAR $

using namespace std;

class Node {

public:
    char cara;
    unordered_map<unsigned int, unsigned int> map_seqId_nbOccs; // hash table: seqId--> list position in this Seq.

    Node* listChildren[NB_CHAR]{nullptr}; //unique_ptr<Node> listChildren[NB_CHAR]{nullptr}
                                            // i think i not need unique_ptr, because i use an object of calss SuffixTree_QuadraritcTime wich is not pointer, at the end (out of scop), the system call delete the object Suff... wich call the destructor, wich call delete Node, wichi is itself recursive call to delete all children :)

    explicit Node(char cara) : cara(cara){}

    // this will delete all the subTree rooted by this node.
    virtual ~Node(){
        for (auto & child : listChildren) {
            delete child;
            child= nullptr;
        }
    }

    //virtual ~Node(){}

    void AddNewPosition(unsigned int idx_seq, unsigned int pos)
    {
        if(this->map_seqId_nbOccs.find(idx_seq) == this->map_seqId_nbOccs.end())
        {
            pair<unsigned int, unsigned int> my_pair = std::make_pair(idx_seq,0);
            my_pair.second=1;

            this->map_seqId_nbOccs.insert(my_pair);
        }
        else
        {
            this->map_seqId_nbOccs[idx_seq]++;
        }
    }

    string toString()
    {
        string str="[";
        str+= this->cara + string(" { ");

        for (auto & child : listChildren)
        {
            if(child)
            {
                str+= child->cara + string(", ");
            }
            else
                str+=string("nul, ");
        }

        str+="}";

        str+=" pos {";


        for (const auto &my_pair:map_seqId_nbOccs)
        {
            cout<<"seqId: "<<my_pair.first;
            cout<<util::to_string(my_pair.second);
        }

        str+="}";

        str+="]";

        return str;
    }

    friend ostream &operator<<(ostream &os, const Node &node)
    {
        os << "[";
        os << node.cara;
        os << "| ";
        for (const auto &my_pair:node.map_seqId_nbOccs)
        {
            os<<"seqId: "<<my_pair.first;
            os<<util::to_string(my_pair.second);
            os<<", ";
        }

        os<<"]";

        return os;
    }
};

#endif //SUFFIX_TREE_NODE_H
