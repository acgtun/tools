/*
 * transfer SAM file to MappedRead format for bsmap
 */

#include <map>
#include <set>
#include <cstdint>
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

void ReadBSMAPResult(const string& file_name,
                     map<string, uint32_t>& read_mismatch) {
  cout << "read bsmap results..." << endl;
  string line;
  string chrom;
  string start_pos;
  string end_pos;
  string read_name;
  int flag;
  int num_of_mismatches;
  string NM;
  string ZS;
  string nothing;
  string read_seq;
  string read_score;

  set<string> reads;

  ifstream fin(file_name.c_str());
  while (getline(fin, line)) {
    if (line[0] == '@')
      continue;
    istringstream iss(line);
    iss >> read_name >> flag >> chrom >> start_pos >> flag >> nothing >> nothing
        >> nothing >> nothing >> read_seq >> read_score >> NM >> ZS;
    if (reads.find(read_name) != reads.end()) {
      read_mismatch[read_name] = 100;
      cout << "multiple ..";
    } else {
      reads.insert(read_name);
      sscanf(NM.c_str(), "NM:i:%d", &num_of_mismatches);
      read_mismatch.insert(make_pair(read_name, num_of_mismatches));
    }
  }
  cout << endl;
}

void ReadResult(const string& file_name, const string& output_file,
                const map<string, uint32_t>& read_mismatch) {
  cout << "transfer..." << endl;
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
  ofstream fout(output_file.c_str());
  while (getline(fin, line)) {
    istringstream iss(line);
    iss >> chrom >> start_pos >> end_pos >> read_name >> num_of_mismatches
        >> strand >> read_seq >> read_score;
    map<string, uint32_t>::const_iterator it = read_mismatch.find(read_name);
    if (it != read_mismatch.end() && it->second != 100) {
      CMAPPINGResult res(chrom, start_pos, end_pos, read_name, it->second,
                         strand, read_seq, read_score);
      res.Output(fout);
    }
  }
  fin.close();
  fout.close();
}

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);

  string bsmap_sam_file, to_mr_file, output_file;
  Option::GetOption("-f", bsmap_sam_file);
  Option::GetOption("-m", to_mr_file);
  Option::GetOption("-o", output_file);

  map<string, uint32_t> read_mismatch;
  ReadBSMAPResult(bsmap_sam_file, read_mismatch);
  ReadResult(to_mr_file, output_file, read_mismatch);

  return 0;
}
