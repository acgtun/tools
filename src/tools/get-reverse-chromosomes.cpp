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

void ReadChromosomes(const string& chrom_file, vector<string>& chrom_names,
                     vector<string>& chrom_seqs) {
  cerr << "[READING CHROMOSOMES] " << endl;

  /* read chromosome name and seqeunce from the chromosome file */
  read_fasta_file(chrom_file.c_str(), chrom_names, chrom_seqs);
}

void ReverseGenome(const vector<string>& chrom_names,
                   const vector<string>& chrom_seqs,
                   const string& output_file) {
  ofstream fout(output_file.c_str());
  for (uint32_t i = 0; i < chrom_seqs.size(); ++i) {
    string reverse_string;
    for (uint32_t j = 0; j < chrom_seqs[i].size(); ++j) {
      reverse_string += complimentBase(
          chrom_seqs[i][chrom_seqs[i].size() - j - 1]);
    }
    cout << chrom_names[i] << endl;
    for (uint32_t j = 0; j < reverse_string.size(); ++j) {
      if (j % 80 == 0 && i != 0)
        fout << endl;
      fout << reverse_string[j] << endl;
    }
  }

  fout.close();
}

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);
  string chrom_file, output_file;
  Option::GetOption("-c", chrom_file);
  Option::GetOption("-o", output_file);

  vector<string> chrom_names;
  vector<string> chrom_seqs;

  ReadChromosomes(chrom_file, chrom_names, chrom_seqs);

  ReverseGenome(chrom_names, chrom_seqs, output_file);

  return 0;
}
