//
// Created by ibra on 5/11/2020.
//

#ifndef MOTIFSEXTRACTIONSELECTION_HASHER_CONTAINER_H
#define MOTIFSEXTRACTIONSELECTION_HASHER_CONTAINER_H

#include <vector>

class hasher_Container{
public:

    // https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
    // boost::hash_combine
    // https://www.boost.org/doc/libs/1_35_0/doc/html/boost/hash_combine_id241013.html

    static unsigned int hash(std::vector<unsigned int> const& vec)
    {
        unsigned int seed = vec.size();
        for(auto& i : vec) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    static unsigned int hash(std::unordered_map<unsigned int, unsigned int> const& umap)
    {
        unsigned int seed = umap.size()*2; // *2 each elemet is a pair: key, value
        for(auto& pair : umap)
        {
            seed ^= pair.first + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= pair.second + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    static unsigned int hash(std::map<unsigned int, unsigned int> const& map)
    {
        unsigned int seed = map.size()*2; // *2 each elemet is a pair: key, value
        for(auto& pair : map)
        {
            seed ^= pair.first + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= pair.second + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    static unsigned long long int hash_64bits(std::vector<unsigned int> const& vec)
    {
        unsigned long long int seed = vec.size();
        for(auto& i : vec) {
            seed ^= i + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    static unsigned long long int hash_64bits(std::unordered_map<unsigned int, unsigned int> const& umap)
    {
        unsigned long long int seed = umap.size()*2; // *2 each elemet is a pair: key, value
        for(auto& pair : umap)
        {
            seed ^= pair.first + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
            seed ^= pair.second + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    static unsigned long long int hash_64bits(std::map<unsigned int, unsigned int> const& map)
    {
        unsigned long long int seed = map.size()*2; // *2 each elemet is a pair: key, value
        for(auto& pair : map)
        {
            seed ^= pair.first + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
            seed ^= pair.second + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
#endif //MOTIFSEXTRACTIONSELECTION_HASHER_CONTAINER_H
