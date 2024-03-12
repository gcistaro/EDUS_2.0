#include <iostream>
#include <fstream>
#include <sstream>
#include <string> 
#include <vector>
/*
 *   ReadFile
 *   Short description:
 *      Read the text file with name FileName splitting it into lines, each line is splitted in words divided by " "
 *   Input:
 *      string& FileName -> Name of the file to read
 *   Output:
 *      file_content[i][j] -> i-th line of file, j-th word
 *
*/
auto ReadFile(const std::string& FileName) -> std::vector<std::vector<std::string>>
{    
    std::ifstream FileName_Stream;
    FileName_Stream.open(FileName.c_str());
    //get all the lines and put them in a container called file_content
    
    std::vector<std::vector<std::string>> file_content;
    std::string line;
    while( std::getline(FileName_Stream,line) ){
        file_content.push_back(std::vector<std::string>());
        std::string aux_word;
        std::stringstream line_stringstream(line);
        while( line_stringstream >> aux_word ){
            file_content[file_content.size()-1].push_back(aux_word);
        }
    }
    FileName_Stream.close();
    return file_content;
}


void WriteFile(const std::string& FileName, const auto& file_content)
{
    std::ofstream FileName_Stream;
    FileName_Stream.open(FileName.c_str());
    
    for(auto& line : file_content){
        for(auto& word : line){
            FileName_Stream << word << " ";
        }
        FileName_Stream << std::endl;
    }
    FileName_Stream.close();

}