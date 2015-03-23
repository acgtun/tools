/*
 * trim-reads.cpp
 *
 *  delete the end of each read, make them shorter
 *
 */

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
#include "smithlab_os.hpp"

using namespace std;

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);

  string fastq_file;
  string trim_fastq_file;
  uint32_t after_read_len;

  Option::GetOption("-f", fastq_file);
  Option::GetOption("-o", trim_fastq_file);
  Option::GetOption("-l", after_read_len, 90);

  vector<string> names;
  vector<string> sequences;
  vector<string> scores;
  read_fastq_file(fastq_file.c_str(), names, sequences, scores);

  ofstream fout(trim_fastq_file.c_str());
  for (uint32_t i = 0; i < names.size(); ++i) {
    fout << "@" << names[i] << endl;
    if (sequences[i].size() < after_read_len) {
      std::cout << "read length: " << sequences[i].size() << endl;
      exit(EXIT_FAILURE);
    }
    fout << sequences[i].substr(0, after_read_len) << endl;
    fout << "+" << names[i] << endl;
    fout << scores[i].substr(0, after_read_len) << endl;
  }
  fout.close();

  return 0;
}

