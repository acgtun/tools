#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "smithlab_os.hpp"

using namespace std;

int main(int argc, const char *argv[]) {
  string files;
  files = argv[1];
  vector<string> file_names;
  if(isdir(files.c_str())) {
    read_dir(files, file_names);
  } else {
    file_names.push_back(files);
  }

  for (uint32_t i = 0; i < file_names.size(); ++i) {
    ifstream fin(file_names[i].c_str());
    string line;
    int line_count = 0;
    bool is_pbs_output_file = false;
    while (getline(fin, line)) {
      line_count++;

      if(line_count == 1)
        continue;
      if(line_count == 2) {
        if(line.find("Begin PBS Prologue") != string::npos) {
          is_pbs_output_file = true;
          continue;
        } else {
          break;
        }
      }

      if (line.find("Resources:") == string::npos)
        continue;

      if(is_pbs_output_file == false) 
        break;
      uint32_t pos = line.find("cput");
      line = line.substr(pos);
      int hour, min, sec;
      long long mem;
      sscanf(line.c_str(), "cput=%d:%d:%d,mem=%lld", &hour, &min, &sec,
             &mem);
      uint32_t num_of_seconds = 0;
      num_of_seconds = hour * 3600 + min * 60 + sec;
      uint32_t p = file_names[i].find_last_of('/');
      if(p != string::npos) {
        file_names[i] = file_names[i].substr(p + 1);
      }
      cout << (double) num_of_seconds / 3600.00 << "\t" 
          << (double) mem / (1024.00 * 1024.00) << "\t"
          <<  file_names[i] << endl; 
    }
    fin.close();
  }
  return 0;
}
