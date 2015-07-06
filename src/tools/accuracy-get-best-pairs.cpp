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
  Option::GetOption("-a2gf", A2Gforward_file);
  Option::GetOption("-a2gr", A2Greverse_file);

  Option::GetOption("-c", genome_file);
  Option::GetOption("-l", read_length, 90);
  Option::GetOption("-l", max_mismatches, 6);
  Option::GetOption("-L", frag_range, 1000);

  map<string, uint32_t> chrom_length;
  ReadGenome(genome_file, chrom_length);

  vector<vector<CMAPPINGResult> > res_1(1000005);
  vector<vector<CMAPPINGResult> > res_2(1000005);

  ReadGroundTruthResults(C2Tforward_file, chrom_length, read_length, res_1);
  ReadGroundTruthResults(C2Treverse_file, chrom_length, read_length, res_1);
  ReadGroundTruthResults(A2Gforward_file, chrom_length, read_length, res_2);
  ReadGroundTruthResults(A2Greverse_file, chrom_length, read_length, res_2);

  for (uint32_t i = 1; i <= 1000000; ++i) {
    sort(res_1[i].begin(), res_1[i].end(), sortCMP);
    sort(res_2[i].begin(), res_2[i].end(), sortCMP);
  }

  /////////////////////////////////////////////
  // FIND THE BEST PAIR
  cout << "find good pair..." << endl;
  vector<CMAPPINGResult> best_result(1000005);
  uint32_t num_of_paried = 0;
  for (uint32_t r = 1; r <= 1000000; ++r) {
    if (r % 10000 == 0) {
      cout << r << endl;
    }
    pair<int, int> best_pair(-1, -1);
    uint32_t min_num_of_mismatch = max_mismatches;
    uint32_t best_times = 0;
    for (uint32_t i = 0; i < res_1[r].size(); ++i) {
      for (uint32_t j = 0; j < res_2[r].size(); ++j) {
        const CMAPPINGResult& r1 = res_1[r][i];
        const CMAPPINGResult& r2 = res_2[r][j];
        if (r1.strand == r2.strand)
          continue;

        uint32_t num_of_mismatch = r1.mismatch + r2.mismatch;
        if (num_of_mismatch > min_num_of_mismatch)
          break;

        if (r1.chrom != r2.chrom) {
          continue;
        }

        uint32_t s1 = r1.start_pos;
        uint32_t s2 = r2.start_pos;
        uint32_t e1 = r1.end_pos;
        uint32_t e2 = r2.end_pos;
        int frag_size = r1.strand == '+' ? (e2 - s1) : (e1 - s2);

        if (frag_size <= 0 || frag_size > frag_range)
          continue;
        if (num_of_mismatch < min_num_of_mismatch) {
          best_pair = make_pair(i, j);
          best_times = 1;
          min_num_of_mismatch = num_of_mismatch;
        } else if (num_of_mismatch == min_num_of_mismatch) {
          best_pair = make_pair(i, j);
          best_times++;
        }
      }
    }

    if (best_times == 1) {
      num_of_paried++;
      const CMAPPINGResult& r1 = res_1[r][best_pair.first];
      best_result[r] = CMAPPINGResult(r1.chrom, r1.start_pos, r1.end_pos,
                                      "QNAME", min_num_of_mismatch, r1.strand,
                                      "SEQ", "SCORE");
    } else {
      best_result[r] = CMAPPINGResult("XXX", 0, 0, "XXX", 100, '#', "XXX",
                                      "XXX");
    }
  }

  ofstream fout("SRR1171450_best.txt");
  for (uint32_t i = 1; i <= 1000000; ++i) {
    const CMAPPINGResult& r = best_result[i];
    fout << r.chrom << "\t" << r.start_pos << "\t" << r.end_pos << "\t"
         << r.read_name << "\t" << r.mismatch << "\t" << r.strand << "\t"
         << r.read_seq << "\t" << r.read_score << endl;

  }
  fout.close();

  return 0;
}
