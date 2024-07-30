#ifndef STREAMFILE_HPP
#define STREAMFILE_HPP

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
auto ReadFile(const std::string& FileName) -> std::vector<std::vector<std::string>>;
void WriteFile(const std::string& FileName, const std::vector<std::vector<std::string>>& file_content);

#endif