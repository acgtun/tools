#include <map>
#include <set>
#include <stdint.h>
#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "smithlab_os.hpp"

using namespace std;

struct QSUB_ERROR {
  QSUB_ERROR(const double& _total_mapping_time,
             const double& _binary_search_time,
             const double& _num_of_full_check, const double& _fullcheck_time,
             const string& _file)
      : total_mapping_time(_total_mapping_time),
        binary_search_time(_binary_search_time),
        num_of_full_check(_num_of_full_check),
        fullcheck_time(_fullcheck_time),
        file(_file) {
  }
  static bool QSUB_OUT_CMP(const QSUB_ERROR& a, const QSUB_ERROR& b) {
    return a.file < b.file;
  }

  void Output() {
    cout << total_mapping_time << "\t" << binary_search_time << "\t"
         << num_of_full_check << "\t" << fullcheck_time << "\t"
         << fullcheck_time / total_mapping_time << "\t" << file << endl;
  }

  double total_mapping_time;
  double binary_search_time;
  double num_of_full_check;
  double fullcheck_time;
  string file;
};

int main(int argc, const char *argv[]) {
  string files;
  files = argv[1];
  vector<string> file_names;
  if (isdir(files.c_str())) {
    read_dir(files, file_names);
  } else {
    file_names.push_back(files);
  }

//  for(uint32_t i = 0;i < file_names.size();++i) {
//    cout << file_names[i] << endl;
//  }
//  cout << "----------------" << endl;

  vector<QSUB_ERROR> qsub_error;
  for (uint32_t i = 0; i < file_names.size(); ++i) {
    if (file_names[i].find("pbs.e") == string::npos)
      continue;
    cout << file_names[i] << endl;
    ifstream fin(file_names[i].c_str());
    string line;
    double mapping_time;
    double binary_search_time;
    double num_of_full_check;
    double fullcheck_time;
    while (getline(fin, line)) {

      uint32_t pos = line.find("[MAPPING TAKES");
      if (pos != string::npos) {
        sscanf(line.c_str(), "[MAPPING TAKES %lf SECONDS]", &mapping_time);
      }

      pos = line.find("GETREGIONTIME:");
      if (pos != string::npos) {
        sscanf(line.c_str(), "GETREGIONTIME: %lf", &binary_search_time);
      }

      pos = line.find("NUMOFFULLCHECKT:");
      if (pos != string::npos) {
        sscanf(line.c_str(), "NUMOFFULLCHECKT: %lf", &num_of_full_check);
      }

      pos = line.find("FULLCHECKTIME:");
      if (pos != string::npos) {
        sscanf(line.c_str(), "FULLCHECKTIME: %lf", &fullcheck_time);
      }

      uint32_t p = file_names[i].find_last_of('/');
      if (p != string::npos) {
        file_names[i] = file_names[i].substr(p + 1);
      }
    }
    qsub_error.push_back(
        QSUB_ERROR((double) mapping_time / 3600.00,
                   binary_search_time / 3600.00, num_of_full_check,
                   fullcheck_time / 3600, file_names[i]));
    fin.close();
  }

  sort(qsub_error.begin(), qsub_error.end(), QSUB_ERROR::QSUB_OUT_CMP);
  for (uint32_t i = 0; i < qsub_error.size(); ++i) {
    qsub_error[i].Output();
  }

  return 0;
}
