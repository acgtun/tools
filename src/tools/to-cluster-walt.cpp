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

void single_end_run(const string& fq_file) {

}

void paired_end_run(const string& fq_file1, const string& fq_file2) {

}

bool check_paired_end(const string& fq_file1, const string& fq_file2) {
  string p1 = fq_file1.find("_1.fastq");
  string p2 = fq_file2.find("_2.fastq");
  if (p1 != string::npos && p2 != string::npos
      && fq_file1.substr(0, p1) == fq_file2.substr(0, p2))
    return true;

  return false;
}

int main(int argc, const char *argv[]) {
  string files;
  files = argv[1];
  vector<string> file_names;
  if (isdir(files.c_str())) {
    read_dir(files, file_names);
  } else {
    file_names.push_back(files);
  }
  sort(file_names.begin(), file_names.end());

  for (uint32_t i = 0; i < file_names.size(); ++i) {
    if (i + 1 >= fele_names.size()) {
      single_end_run(file_names[i]);
    } else {

    }

  }

  return 0;
}
