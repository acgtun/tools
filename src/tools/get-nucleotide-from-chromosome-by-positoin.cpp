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


/* get the compliment strand nucleotide */
inline char complimentBase(const char& nt) {
  switch (nt) {
    case 'a':
      return ('t');
    case 'c':
      return ('g');
    case 'g':
      return ('c');
    case 't':
      return ('a');
    case 'A':
      return ('T');
    case 'C':
      return ('G');
    case 'G':
      return ('C');
    case 'T':
      return ('A');
    default:
      return ('N');
  }
}


int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);

  string chrom_file;
  uint32_t chrom_pos;
  uint32_t read_len;
  bool rc = false;
  Option::GetOption("-c", chrom_file);
  Option::GetOption("-p", chrom_pos, 0);
  Option::GetOption("-l", read_len, 100);
  Option::ChkStrExist("-r", rc);


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

  for(uint32_t i = 0;i < chrom_seq.size();++i) {
    chrom_seq[i] = toupper(chrom_seq[i]);
  }

  if(rc) {
    string chrom_seq_rc;
    for (uint32_t j = 0; j < chrom_seq.size(); ++j) {
      chrom_seq_rc += complimentBase(
          chrom_seq[chrom_seq.size() - j - 1]);
    }
    cout << "size of chromosome " << chrom_file << " is " << chrom_seq.size() << endl;
    cout << "reverse strand" << endl;
    cout << "start position: " << chrom_pos << endl;
    uint32_t i, j;
    uint32_t start = chrom_seq.size() - chrom_pos - read_len;
    for (i = start, j = 0; j < 90; ++i, ++j) {
      cout << chrom_seq_rc[i];
    }
    cout << endl;
  } else {
    cout << "size of chromosome " << chrom_file << " is " << chrom_seq.size() << endl;
    cout << "start position: " << chrom_pos << endl;
    for (uint32_t i = chrom_pos, j = 0; i < chrom_seq.size() && j < read_len; ++i, ++j) {
      cout << chrom_seq[i];
    }
    cout << endl;
  }

  return 0;
}
