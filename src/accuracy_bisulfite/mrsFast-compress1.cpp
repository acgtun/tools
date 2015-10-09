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

bool IsForwardStrand(const uint32_t& flag) {
  return !(flag & 0x10);
}

void SAMfileReader(const string& result_file, const string& strand) {
  ifstream fin(result_file.c_str());
  string line;
  string qname, rname, cigar, rnext, seq, qual;
  uint32_t flag, pos, mapq, pnext, tlen;
  string NM, MD;
  ofstream fout(string(result_file + "_post_processing.txt"));
  while (getline(fin, line)) {
    if (line[0] == '@')
      continue;
    istringstream iss(line);
    iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext
        >> tlen >> seq >> qual >> NM >> MD;
    if (!IsForwardStrand(flag)) {
      continue;
    }
    int mismatch;
    sscanf(NM.c_str(), "NM:i:%d", &mismatch);
    fout << qname << " " << rname << " " << pos << " " << strand << " "
         << mismatch << endl;
  }
  fin.close();
  fout.close();
}

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);
  string result_file;
  int st;
  string strand;
  Option::GetOption("-f", result_file);
  Option::GetOption("-s", st, 3);

  strand = st == 0 ? '+' : '-';
  SAMfileReader(result_file, strand);

  return 0;
}
