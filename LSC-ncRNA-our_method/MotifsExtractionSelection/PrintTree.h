//
// Created by ibra on 11/29/2019.
//

#ifndef SUFFIX_TREE_PRINTTREE_H
#define SUFFIX_TREE_PRINTTREE_H

#include <string>
#include <iostream>

#include "Node.h"

using namespace std;


class PrintTree {

    using bits=std::vector<bool>;

    static void p_tabs(const bits & b)
    {
        for (auto x: b)
            //std::cout << (x?" \u2502":"  ");
            std::cout << (x?" |":"  ");
    }

    static void p_show(const Node *r, bits &b)
    {
        // https://en.wikipedia.org/wiki/Box-drawing_character
        if (r)
        {
            std::cout
                    //<< "[" << r->cara << "," << util::to_string(*r->listPosition) << "]"
                    //<< "[" << r->toString() << "]"
                    << *r
                    << std::endl;


            for (int i = 0; i <NB_CHAR-1 ; ++i)
            {
                p_tabs(b);
                ///std::cout << " \u251c"; // ├
                std::cout << " |--"; // ├
                b.push_back(true);
                p_show(r->listChildren[i],b);
                b.pop_back();
            }

            p_tabs(b);
            //std::cout << " \u2514"; // └
            std::cout << " L-"; // └
            b.push_back(false);
            p_show(r->listChildren[NB_CHAR-1],b);
            b.pop_back();
        }
        else
            //std::cout << " \u25cb" << std::endl; // ○
            std::cout << " o" << std::endl; // ○
    }

    static void show(const Node *root)
    {
        bits b;
        p_show(root,b);
    }

public:

    static void PrintSuffixTree(const Node *rootTree)
    {
        show(rootTree);
    }

};


#endif //SUFFIX_TREE_PRINTTREE_H
