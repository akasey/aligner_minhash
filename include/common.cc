//
// Created by Akash Shrestha on 5/7/18.
//

#include "common.h"

void write_in_file(const void *ptr, size_t size, size_t count, FILE *stream) {
    if (fwrite(ptr, size, count, stream) != count) {
        LOG(ERROR) << "Error writing into file";
        exit(-10);
    }
}

void read_from_file(void *ptr, size_t size, size_t count, FILE *stream) {
    if (fread(ptr, size, count, stream) != count) {
        LOG(ERROR) << "Error reading from file";
        exit(-11);
    }
}

std::string & ltrim(std::string & str) {
    auto it2 =  std::find_if( str.begin() , str.end() , [](char ch){ return !std::isspace<char>(ch , std::locale::classic() ) ; } );
    str.erase( str.begin() , it2);
    return str;
}

std::string & rtrim(std::string & str) {
    auto it1 =  std::find_if( str.rbegin() , str.rend() , [](char ch){ return !std::isspace<char>(ch , std::locale::classic() ) ; } );
    str.erase( it1.base() , str.end() );
    return str;
}

std::string & trim(std::string & str) {
    return ltrim(rtrim(str));
}

bool str_endswith(std::string &original, std::string &with) {
    std::size_t found = original.rfind(with);
    if (found==std::string::npos || found != original.length()-with.length())
        return false;
    return true;
}

void mkdir(std::string dirname) {
    const int dir_err = system(("mkdir -p "+dirname).c_str());
    if (-1 == dir_err)
    {
        LOG(ERROR) << "Couldn't mkdir " + dirname;
        exit(-1);
    }
}

std::vector<std::string> getFilesInDirectory(std::string directory, std::string extensionFilter) {
    DIR *dirp = opendir(directory.c_str());
    if(!dirp) {
        LOG(ERROR) << "Could open directory " << directory;
        exit(-33);
    }
    std::vector<std::string> toReturn;
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        std::string filename = std::string(dp->d_name);
        if (filename.compare(".")==0 || filename.compare("..")==0 || (!extensionFilter.empty() && !str_endswith(filename, extensionFilter)))
            continue;
        std::string fqdn = directory + "/" + std::string(dp->d_name);
        toReturn.push_back(fqdn);
    }
    return toReturn;
}

std::vector<std::string> split(std::string sentence, std::string delimiter){
    std::vector<std::string> list;
    size_t pos = 0;
    std::string token;
    while ((pos = sentence.find(delimiter)) != std::string::npos) {
        token = sentence.substr(0, pos);
        list.push_back(token);
        sentence.erase(0, pos + delimiter.length());
    }
    list.push_back(sentence);
    return list;
}

std::string basename(std::string &fqdn) {
    std::size_t found = fqdn.rfind("/");
    if (found==std::string::npos)
        found = -1;
    return fqdn.substr(found+1);
}

int absoluteValue(int a) {
    return a<0 ? -a : a;
}

float absoluteValue(float a) {
    return a<0 ? -a : a;
}