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

struct CMAPPINGResult {
  CMAPPINGResult(string _chrom, string _start_pos, string _end_pos,
                 string _read_name, int _num_of_mismatches, char _strand,
                 string _read_seq, string _read_score) {
    chrom = _chrom;
    start_pos = _start_pos;
    end_pos = _end_pos;
    read_name = _read_name;
    num_of_mismatches = _num_of_mismatches;
    strand = _strand;
    read_seq = _read_seq;
    read_score = _read_score;
  }

  void Output(ofstream& fout) {
    fout << chrom << "\t" << start_pos << "\t" << end_pos << "\t" << read_name
         << "\t" << num_of_mismatches << "\t" << strand << "\t" << read_seq
         << "\t" << read_score << endl;
  }

  string chrom;
  string start_pos;
  string end_pos;
  string read_name;
  int num_of_mismatches;
  char strand;
  string read_seq;
  string read_score;
};

void ReadResult(const string& file_name, map<string, CMAPPINGResult>& res) {
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
    res.insert(
        make_pair(
            read_name,
            CMAPPINGResult(chrom, start_pos, end_pos, read_name,
                           num_of_mismatches, strand, read_seq, read_score)));
  }
}

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);

  string mapping_file, output_file;
  Option::GetOption("-f", mapping_file);
  Option::GetOption("-o", output_file);

  map<string, CMAPPINGResult> res;
  ReadResult(mapping_file, res);

  map<uint32_t, uint32_t> mismatch_count;
  for (map<string, CMAPPINGResult>::const_iterator it = res.begin();
      it != res.end(); ++it) {
    mismatch_count[it->second.num_of_mismatches]++;
  }

  ofstream fout(output_file.c_str());
  for (map<uint32_t, uint32_t>::const_iterator it = mismatch_count.begin();
      it != mismatch_count.end(); ++it) {
    fout << it->first << " " << it->second << endl;
  }
  fout.close();

  return 0;
}

