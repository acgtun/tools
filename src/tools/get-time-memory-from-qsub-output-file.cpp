#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "option.hpp"

using namespace std;

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);
  int num_of_files;
  Option::GetOption("-n", num_of_files, 1);

  vector<string> files(num_of_files);
  for (int i = 0; i < num_of_files; ++i) {
    char file[4];
    sprintf(file, "-f%d", i + 1);
    Option::GetOption(file, files[i]);
  }

  for (int i = 0; i < num_of_files; ++i) {
    ifstream fin(files[i].c_str());
    string line;
    while (getline(fin, line)) {
      if (line.find("Resources:") == string::npos)
        continue;
      uint32_t pos = line.find("cput");
      line = line.substr(pos);
      int hour, min, sec;
      long long mem;
      sscanf(line.c_str(), "cput=%d:%d:%d,mem=%lld", &hour, &min, &sec,
             &mem);
      uint32_t num_of_seconds = 0;
      num_of_seconds = hour * 3600 + min * 60 + sec;
      uint32_t p = files[i].find_last_of('/');
      if(p != string::npos) {
        files[i] = files[i].substr(p + 1);
      }
      cout << (double) num_of_seconds / 3600.00 << "\t" 
          << (double) mem / (1024.00 * 1024.00) << "\t"
          <<  files[i] << endl; 
    }
  }
  return 0;
}
