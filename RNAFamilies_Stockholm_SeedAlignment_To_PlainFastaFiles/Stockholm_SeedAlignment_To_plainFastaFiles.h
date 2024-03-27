//
// Created by ibra on 2/11/2020.
//

#ifndef RNAFAMILIES_STOCKHOLM_SEEDALIGNMENT_TO_PLAINFASTAFILES_STOCKHOLM_SEEDALIGNMENT_TO_PLAINFASTAFILES_H
#define RNAFAMILIES_STOCKHOLM_SEEDALIGNMENT_TO_PLAINFASTAFILES_STOCKHOLM_SEEDALIGNMENT_TO_PLAINFASTAFILES_H

#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>


class Stockholm_SeedAlignment_To_plainFastaFiles{

    static const char GAP='-';

    // https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
    // trim from start
    static inline std::string &ltrim(std::string &s)
    {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                        std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    // trim from end
    static inline std::string &rtrim(std::string &s)
    {
        s.erase(std::find_if(s.rbegin(), s.rend(),
                             std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    // trim from both ends
    static inline std::string &trim(std::string &s)
    {
        return ltrim(rtrim(s));
    }

    static std::string getUnGapedStr(const std::string &str)
    {
        std::string ungaped_str;

        for (char c : str)
        {
            if(c!=GAP) ungaped_str+=c;
        }

        return ungaped_str;
    }

    static bool saveToFile(const std::string & file_name, const std::string &buffer)
    {
        std::ofstream output;
        output.open (file_name);

        if(!output.good())
            return false;

        output<<buffer;

        output.close();

        return true;
    }

public:

    static uint32_t convert(const std::string &file_name, const std::string &dir_output)
    {
        std::ifstream input(file_name);

        if(!input.good()){
            std::cout<<"can't read the file, input is not good ..."<<std::endl;
            return 0; // return 0 family are treated
        }
            

        std::string line;
        std::string name; // > id of seq
        std::string content; // seq

        std::string buffer;

        uint32_t nb_families=0;

        const std::string stk_code_rfam_id = "#=GF AC";
        const std::string stk_code_RNA_id = "#=GF ID";
        const std::string stk_code_end_alignment = "//";
        std::string rfam_id;
        std::string RNA_id;
        std::string RNA_file_name;


        std::cout<<" file opened ok "<<std::endl;

        while( std::getline( input, line ).good() )
        {
            if(line.empty()) continue;

            if(line[0] == '#' ) // Identifier marker for a new entry
            {
                //std::cout<<"in ## line 0 = "<<line[0]<<std::endl;

                if(line.rfind(stk_code_rfam_id,0)==0)
                {
                    //std::cout<<"stk_code_rfam_id.size() : "<<stk_code_rfam_id.size()<<std::endl;

                    rfam_id = line.substr(stk_code_rfam_id.size());
                    trim(rfam_id);
                    nb_families++;
                }
                else if(line.rfind(stk_code_RNA_id,0)==0)
                {
                    //std::cout<<"stk_code_RNA_id.size() : "<<stk_code_RNA_id.size()<<std::endl;
                    RNA_id = line.substr(stk_code_RNA_id.size());
                    trim(RNA_id);

                    RNA_file_name+=dir_output;
                    RNA_file_name+="/";
                    RNA_file_name +=rfam_id;
                    /// RNA_file_name +="__";
                    /// RNA_file_name +=RNA_id;  // can be add (the RNA family name) it is a choice ...
                    RNA_file_name +=".fasta";
                    RNA_file_name +=".txt";   // can be removed , it a choice ...

                    //std::cout<<"name : "<<RNA_file_name<<std::endl;
                }
            }
            else if(line[0]=='/')//  && line[1]=='/';  // end of alignment
            {
                std::cout<<" save File :"<<RNA_file_name<<std::endl;
                saveToFile(RNA_file_name,buffer);
                buffer.clear();
                RNA_file_name.clear();
            }
            else if(std::isspace(line[0])) // line[0]=='\n' , \r , \r\n
            {
                continue;
            }
            else
            {
                char space=' ';
                uint32_t pos_space = line.find(space);

                name = line.substr(0,pos_space);
                uint32_t pos_name_separator = name.find('/');
                name[pos_name_separator] = '_';

                content = line.substr(pos_space);
                trim(content);
                content=getUnGapedStr(content);

                buffer+='>';
                buffer+=name;
                buffer+="\n";
                buffer+=content;
                buffer+="\n";
            }
        }

        input.close();

        return nb_families;
    }
};
#endif //RNAFAMILIES_STOCKHOLM_SEEDALIGNMENT_TO_PLAINFASTAFILES_STOCKHOLM_SEEDALIGNMENT_TO_PLAINFASTAFILES_H
