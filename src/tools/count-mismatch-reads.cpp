/*
 * count number of matched reads for different number of mismatches
 * unordered_map<uint32_t, uint32_t> mismatch_reads;
 * key: number of mismatches
 * value: number of matched reads
 */
#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "option.hpp"

using namespace std;

void MismatchCount(const string& file_name,
                   map<uint32_t, uint32_t>& mismatch_count) {
  string line;
  string chrom;
  string start_pos;
  string end_pos;
  string read_name;
  int num_of_mismatches;
  char strand;
  string read_seq;
  string read_score;

  ifstream fin(file_name.c_str());
  while (getline(fin, line)) {
    istringstream iss(line);
    iss >> chrom >> start_pos >> end_pos >> read_name >> num_of_mismatches
        >> strand >> read_seq >> read_score;
    mismatch_count[num_of_mismatches]++;
  }
}

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);

  string mapping_file, output_file;
  Option::GetOption("-f", mapping_file);
  Option::GetOption("-o", output_file);

  map<uint32_t, uint32_t> mismatch_count;
  MismatchCount(mapping_file, mismatch_count);

  ofstream fout(output_file.c_str());
  for (map<uint32_t, uint32_t>::const_iterator it = mismatch_count.begin();
      it != mismatch_count.end(); ++it) {
    fout << it->first << " " << it->second << endl;
  }
  fout.close();

  return 0;
}
