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

  string chrom_file;
  uint32_t chrom_pos;
  Option::GetOption("-c", chrom_file);
  Option::GetOption("-p", chrom_pos, 0);

  ifstream fin(chrom_file.c_str());
  string line;
  string chrom_seq;
  while (getline(fin, line)) {
    if (line[0] == '>')
      continue;
    else {
      chrom_seq += line;
    }
  }
  cout << "size of chromosome " << chrom_file << " is " << chrom_seq.size() << endl;
  cout << "start position: " << chrom_pos << endl;
  for (uint64_t i = chrom_pos, j = 0; i < chrom_seq.size() && j < 100; ++i, ++j) {
    cout << chrom_seq[i];
  }
  cout << endl;
  return 0;
}
