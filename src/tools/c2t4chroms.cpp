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

void ReadChromosomes(const string& chrom_file, vector<string>& chrom_names,
                     vector<string>& chrom_seqs) {
  cerr << "[READING CHROMOSOMES] " << endl;

  /* read chromosome name and seqeunce from the chromosome file */
  read_fasta_file(chrom_file.c_str(), chrom_names, chrom_seqs);
}

char C2T(const char& chr, const bool& AGWILDCARD) {
  if (AGWILDCARD) {
    if (chr == 'G' || chr == 'g')
      return 'A';
    return chr;
  } else {
    if (chr == 'C' || chr == 'c')
      return 'T';
    return chr;
  }
}

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

string ReverseString(const string& str) {
  uint32_t size = str.size();
  string ret(size, 'N');
  for (uint32_t i = 0; i < size; ++i) {
    ret[i] = str[size - i - 1];
  }

  return ret;
}

string ReverseComplimentString(const string& str) {
  string ret = ReverseString(str);
  for (uint32_t i = 0; i < str.size(); ++i) {
    ret[i] = complimentBase(ret[i]);
  }

  return ret;
}

void C2TChromosomes(const vector<string>& chrom_names,
                    const vector<string>& chrom_seqs, const string& output_file,
                    const bool& bnewline, const bool& AGWILDCARD,
                    const bool& reverse_compliment) {
  ofstream fout(output_file.c_str());
  for (uint32_t i = 0; i < chrom_seqs.size(); ++i) {
    string reverse_string;
    string chrom_string = chrom_seqs[i];
    if (reverse_compliment) {
      chrom_string = ReverseComplimentString(chrom_string);
    }
    for (uint32_t j = 0; j < chrom_string.size(); ++j) {
      if (j % 50 == 0 && j != 0 && bnewline) {
        reverse_string += '\n';
      }
      //cout << "(" << chrom_string[j] << "," << C2T(chrom_string[j], AGWILDCARD) << ")" << endl;
      reverse_string += C2T(chrom_string[j], AGWILDCARD);
    }
    fout << ">" << chrom_names[i] << endl;
    cout << chrom_names[i] << endl;
    if (reverse_string.size() != 0
        && reverse_string[reverse_string.size() - 1] != '\n') {
      reverse_string += '\n';
    }

    fout << reverse_string;
  }

  fout.close();
}

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);
  bool bnewline = false;
  bool AGWILDCARD = false;
  bool reverse_compliment = false;
  string chrom_file, output_file;
  Option::GetOption("-c", chrom_file);
  Option::GetOption("-o", output_file);
  Option::ChkStrExist("-A", AGWILDCARD);
  Option::ChkStrExist("-newline", bnewline);
  Option::ChkStrExist("-rc", reverse_compliment);

  vector<string> chrom_names;
  vector<string> chrom_seqs;

  ReadChromosomes(chrom_file, chrom_names, chrom_seqs);
  //RevChromosomes(chrom_names, chrom_seqs);

  //return 0;

  C2TChromosomes(chrom_names, chrom_seqs, output_file, bnewline, AGWILDCARD,
                 reverse_compliment);

  return 0;
}
