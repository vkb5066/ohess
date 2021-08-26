#pragma once

#include <iostream>
#include <string>
#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST

#if defined(_WIN32) //also includes windows 64
#include <direct.h>   
#define getcwd _getcwd // stupid MSFT "deprecation" warning
#else
#include <unistd.h>
#endif

//Various utilities that are meant to work on both windows and linux

/*
* Delimiter for file navigation
*/
#if defined(_WIN32) 
const extern std::string dl = "\\";
#else
const extern std::string dl = "/";
#endif

/*
* Get current working directory 
* From https://stackoverflow.com/questions/2868680/what-is-a-cross-platform-way-to-get-the-current-directory
*/
std::string GetCwd(){
    return *new std::string(getcwd(NULL, 0));
}


/*
* Create directory path
* From 'https://stackoverflow.com/questions/675039/how-can-i-create-directory-tree-in-c-linux'
*/
bool DoesDirExist(const std::string& path){
#if defined(_WIN32)
    struct _stat info;
    if (_stat(path.c_str(), &info) != 0)
    {
        return false;
    }
    return (info.st_mode & _S_IFDIR) != 0;
#else 
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
    {
        return false;
    }
    return (info.st_mode & S_IFDIR) != 0;
#endif
}

bool MakePath(const std::string& path){
#if defined(_WIN32)
    int ret = _mkdir(path.c_str());
#else
    mode_t mode = 0755;
    int ret = mkdir(path.c_str(), mode);
#endif
    if (ret == 0)
        return true;

    switch (errno){
        case ENOENT:{
        // parent didn't exist, try to create it
            int pos = path.find_last_of('/');
            if (pos == std::string::npos)
#if defined(_WIN32)
                pos = path.find_last_of('\\');
            if (pos == std::string::npos)
#endif
                return false;
            if (!MakePath( path.substr(0, pos) ))
                return false;
        }
        // now, try to create again
#if defined(_WIN32)
        return 0 == _mkdir(path.c_str());
#else 
        return 0 == mkdir(path.c_str(), mode);
#endif

    case EEXIST:
        return DoesDirExist(path);

    default:
        return false;
    }
}