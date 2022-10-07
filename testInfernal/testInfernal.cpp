#include <iostream>
#include <string>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <climits>
#include <unistd.h>
#include <chrono>
#include <iomanip>
#include <fstream>


using namespace std;

inline char separator()
{
#ifdef _WIN32
    return '\\';
#else
    return '/';
#endif
}


std::string getFileName(const std::string& filePath, bool withExtension = true, char seperator = separator())
{
    // Get last dot position
    std::size_t dotPos = filePath.rfind('.');
    std::size_t sepPos = filePath.rfind(seperator);

    if(sepPos != std::string::npos)
    {
        return filePath.substr(sepPos + 1, filePath.size() - (withExtension || dotPos != std::string::npos ? 1 : dotPos) );
    }
    return "";
}

vector<string> getListFiles(const string &path_dir)
{
    vector<string> list_files;

    DIR *dir;
    struct dirent *entry;
    struct stat file_stat{};

    dir = opendir (path_dir.c_str());

    if(dir== nullptr)
    {
        perror ("Could not open directory");
        return list_files; // return an empty list
    }

    string path_tmp;

    // get all the files within directory
    while ((entry = readdir (dir)) != nullptr)
    {
        path_tmp=path_dir+separator()+entry->d_name;

        if( stat(path_tmp.c_str(),&file_stat) <0) // get the stat
        {
            perror("stat error");
            continue; // go to the next entry
        }

        if (S_ISREG(file_stat.st_mode)) {
            //list_files.emplace_back(entry->d_name);
            list_files.emplace_back(path_tmp); // all the path not just the name of the file
        }
    }

    closedir (dir);

    return list_files;
}

void multiple_alignment(string &dir_input, string &dir_output, string &test_name)
{
    auto list_files_paths = getListFiles(dir_input);

    string command;
    string file_name;

    for (auto & file_path : list_files_paths)
    {
        file_name = getFileName(file_path);

        if(file_name == "info_test") continue;

        //command += "printf \"";
  //      command += "./clustalo -i ";
        command += "clustalo -i ";
        command += file_path;
        command += " -o ";
        command += dir_output;
        command += separator();
        command += "al_"+test_name;
        command += "_"+file_name;
        command += " --outfmt=st --threads=1 -v --force";
        command += " > output_"+test_name;
        //command += " \n \"";

        system(command.c_str());
        command.clear();
    }
}

void secondary_structur(string &dir_input, string &test_name)
{
    auto list_files_paths = getListFiles(dir_input);

    string command;
    string file_name;

    // be carful, RNAalifold output in the same working dirctory where it is lenched

    for (auto & file_path : list_files_paths)
    {
        file_name = getFileName(file_path);

        //command += "printf \"";
        command += "RNAalifold --mis -q --noPS ";
        command += file_path;
        command += " --aln-stk=se_"+file_name;
        command += ".sto";
        command += " > output_"+test_name;
        //command += " \n \"";

        system(command.c_str());
        command.clear();
    }
}

//void cmbuild(string &dir_input, string &test_name) // files are in the same dirctory
void cmbuild(string &test_name) // files are in the same dirctory
{

    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != nullptr)
    {
        printf("Current working dir: %s\n", cwd);
    } else {
        perror("getcwd() error");
        return;
    }


    //auto list_files_paths = getListFiles(dir_input);
    auto list_files_paths = getListFiles(cwd);

    string command;
    string file_name;

    for (auto & file_path : list_files_paths)
    {
        file_name = getFileName(file_path);

        if(file_name.rfind("se_al_"+test_name,0)==0) {
            //command += "printf \"";
            command += "cmbuild ";
            command += file_name+".cm ";;
            command += file_path;
            command += " > output_"+test_name;
            //command += " \n \"";

            system(command.c_str());
            command.clear();
        }
    }
}

void combine_cms_in_one_file(string &test_name) // files are in the same dirctory
{
    string command;

    command += "cat se_al_";
    command += test_name;
    command += "_*.cm > cms_"+test_name+".cm";

    system(command.c_str());
}

void delete_temprary_files(string &dir_input, string &dir_output, string &test_name)
{

    //system(string("rm -r "+dir_output).c_str());
}

void cmsearch_cms_in_one_file(string &dir_input, string &dir_output, string &file_output)
{

}

void time_calculation(std::chrono::system_clock::time_point start, std::chrono::system_clock::time_point end, const string &work_name, const string &test_name)
{

    ofstream file_output;
    file_output.open("time_"+test_name, std::ios_base::app);


    // Calculating total time taken by the program.
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();

    time_taken *= 1e-9;

    file_output << "Time taken by program ( " <<work_name<<" ) is = "<< fixed
         << time_taken << setprecision(9);
    file_output << " sec" << endl;

    file_output.close();
}

int main(int argc, char** argv)
{
    string dir_input = argv[1];
    string dir_output = argv[2];
    string test_name = argv[3];

    std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point end;

    start = chrono::high_resolution_clock::now();
    multiple_alignment(dir_input,dir_output,test_name);
    end = chrono::high_resolution_clock::now();
    time_calculation(start,end,"multiple_alignment", test_name);

    start = chrono::high_resolution_clock::now();
    secondary_structur(dir_output,test_name);
    end = chrono::high_resolution_clock::now();
    time_calculation(start,end,"secondary_structur",test_name);


    start = chrono::high_resolution_clock::now();
    cmbuild(test_name);
    end = chrono::high_resolution_clock::now();
    time_calculation(start,end,"cmbuild",test_name);

    combine_cms_in_one_file(test_name);

    delete_temprary_files(dir_input,dir_output,test_name);

    return 0;
}