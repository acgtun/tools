#include <map>
#include <set>
#include <limits>
#include <vector>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "option.hpp"
#include "smithlab_os.hpp"

using namespace std;

const int MAX_LINE_LENGTH = 10000;
const int MAX_INT = std::numeric_limits<int>::max();

struct CMAPPINGResult {
  CMAPPINGResult(string _chrom = "haoma", uint32_t _start_pos = 0,
                 uint32_t _end_pos = 0, string _read_name = "haoma",
                 int _mismatch = 0, char _strand = '+', string _read_seq =
                     "haoma",
                 string _read_score = "haoma") {
    chrom = _chrom;
    start_pos = _start_pos;
    end_pos = _end_pos;
    read_name = _read_name;
    mismatch = _mismatch;
    strand = _strand;
    read_seq = _read_seq;
    read_score = _read_score;
  }
  friend bool CheckDifference(const CMAPPINGResult& r1,
                              const CMAPPINGResult& r2) {
    return r1.chrom == r2.chrom && r1.start_pos == r2.start_pos
        && r1.end_pos == r2.end_pos && r1.read_name == r2.read_name
        && r1.mismatch == r2.mismatch && r1.strand == r2.strand
        && r1.read_seq == r2.read_seq && r1.read_score == r2.read_score;
  }

  void Output(ofstream& fout) {
    fout << chrom << "\t" << start_pos << "\t" << end_pos << "\t" << read_name
         << "\t" << mismatch << "\t" << strand << "\t" << read_seq << "\t"
         << read_score << endl;
  }
  string chrom;
  uint32_t start_pos;
  uint32_t end_pos;
  string read_name;
  int mismatch;
  char strand;
  string read_seq;
  string read_score;
};

void ReadGenome(const string& genome_file,
                map<string, uint32_t>& chrom_length) {
  fprintf(stderr, "[READING CHROMOSOMES]\n");
  vector<string> chrom_names;
  vector<string> chrom_seqs;

  /* read chromosome name and seqeunce from the chromosome file */
  read_fasta_file(genome_file.c_str(), chrom_names, chrom_seqs);
  chrom_length.clear();
  for (uint32_t i = 0; i < chrom_names.size(); ++i) {
    chrom_length.insert(make_pair(chrom_names[i], chrom_seqs[i].size()));
  }
}

void ReadGroundTruthResults(const string& file_name,
                            const map<string, uint32_t>& chrom_length,
                            const uint32_t read_length,
                            vector<vector<CMAPPINGResult> >& res) {
  cout << "read ground truth results... " << file_name << endl;
  FILE * fin = fopen(file_name.c_str(), "r");
  uint32_t read, num_of_pos;
  uint32_t pos, mismatch;
  uint32_t start_pos, end_pos;
  char strand;
  char chrom[MAX_LINE_LENGTH];
  while (fscanf(fin, "%u %u", &read, &num_of_pos) != EOF) {
    //cout << read << num_of_pos << endl;
    //cerr << "gt "<< read << " " << num_of_pos << endl;
    for (uint32_t i = 0; i < num_of_pos; ++i) {
      fscanf(fin, "%s %u %c %u", chrom, &pos, &strand, &mismatch);
      map<string, uint32_t>::const_iterator it = chrom_length.find(chrom);
      start_pos = pos - 1;
      start_pos =
          strand == '+' ? start_pos : it->second - start_pos - read_length;
      end_pos = start_pos + read_length;

      res[read].push_back(
          CMAPPINGResult(chrom, start_pos, end_pos, "QNAME", mismatch, strand,
                         "SEQ", "SCORE"));
    }
  }

  fclose(fin);
}

bool sortCMP(const CMAPPINGResult& a, const CMAPPINGResult& b) {
  return a.mismatch < b.mismatch;
}

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);

  string C2Tforward_file;
  string C2Treverse_file;
  string A2Gforward_file;
  string A2Greverse_file;
  string genome_file;
  uint32_t read_length;
  uint32_t max_mismatches;
  int frag_range;

  string mapper;
  string mapper_result_file;

  Option::GetOption("-c2tf", C2Tforward_file);
  Option::GetOption("-c2tr", C2Treverse_file);

  Option::GetOption("-c", genome_file);
  Option::GetOption("-l", read_length, 90);
  Option::GetOption("-m", max_mismatches, 6);
  Option::GetOption("-L", frag_range, 1000);

  map<string, uint32_t> chrom_length;
  ReadGenome(genome_file, chrom_length);

  vector<vector<CMAPPINGResult> > res_1(1000005);
  vector<vector<CMAPPINGResult> > res_2(1000005);

  ReadGroundTruthResults(C2Tforward_file, chrom_length, read_length, res_1);
  ReadGroundTruthResults(C2Treverse_file, chrom_length, read_length, res_1);

  for (uint32_t i = 1; i <= 1000000; ++i) {
    sort(res_1[i].begin(), res_1[i].end(), sortCMP);
  }

  /////////////////////////////////////////////
  // FIND THE BEST PAIR
  ofstream fout("best_signle_benchmark.txt");
  for (uint32_t r = 1; r <= 1000000; ++r) {
    if (r % 10000 == 0) {
      cout << r << endl;
    }
    if (res_1[r].size() == 0) {
      fout << "XXX" << "\t" << 0 << "\t" << 100 << endl;

    } else if (res_1[r].size() == 1) {
      fout << res_1[r][0].chrom << "\t" << res_1[r][0].start_pos << "\t"
           << res_1[r][0].mismatch << endl;
    } else {
      if (res_1[r][1].mismatch == res_1[r][0].mismatch) {
        fout << "XXX" << "\t" << 0 << "\t" << 100 << endl;
      } else {
        fout << res_1[r][0].chrom << "\t" << res_1[r][0].start_pos << "\t"
             << res_1[r][0].mismatch << endl;
      }
    }
  }

  fout.close();

  return 0;
}
