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

  void Output() const {
    cout << chrom << "\t" << start_pos << "\t" << end_pos << "\t" << read_name
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

void CompareMappingResultsMR(const string& file_name,
                             const vector<CMAPPINGResult>& best_results,
                             const uint32_t& num_of_paried) {
  cerr << "CompareMappingResultsMR..." << endl;
  cerr << file_name << endl;
  string line;
  string chrom;
  uint32_t start_pos;
  uint32_t end_pos;
  string read_name;
  int num_of_mismatches;
  char strand;
  string read_seq;
  string read_score;
  uint32_t cur_num_of_paried = 0;
  unsigned int read;
  cerr << file_name << endl;
  ifstream fin(file_name);
  while (getline(fin, line)) {
    istringstream iss(line);
    iss >> chrom >> start_pos >> end_pos >> read_name >> num_of_mismatches
        >> strand >> read_seq >> read_score;
    if (read_name.substr(0, 5) == "FRAG:") {
      sscanf(read_name.c_str(), "FRAG:SRR1171540.%u", &read);
      if (read > 1000000)
        break;

      if (best_results[read].chrom == "XXX"
          && best_results[read].mismatch == 100) {
        continue;
      }

      if (num_of_mismatches == best_results[read].mismatch
          && chrom == best_results[read].chrom) {
        cur_num_of_paried++;
      } else {
        cout << "FRAG:SRR1171540." << read << endl;
        best_results[read].Output();
        cout << line << endl;
        cout << "---------------------------------------------------" << endl;
      }
    }
  }
  cout << cur_num_of_paried << "\t" << num_of_paried << "\t"
       << static_cast<double>(cur_num_of_paried) / num_of_paried << endl;

  fin.close();
}

void ReadBsmapResults(const char* file_name, vector<CMAPPINGResult>& res) {
  FILE * fin = fopen(file_name, "r");
  if (!fin) {
    cerr << "cannot open input file " << file_name << endl;
  }
  string qname, rname, cigar, rnext, seq, qual, NM, MD;
  uint32_t flag, mapq, pnext, tlen, pos, end_pos;
  int mismatch;
  char strand = '+';
  char cline[MAX_LINE_LENGTH];
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    if (cline[0] == '@')
      continue;
    istringstream iss(cline);
    iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext
        >> tlen >> seq >> qual >> NM >> MD;
    sscanf(NM.c_str(), "NM:i:%d", &mismatch);
    unsigned int read;
    sscanf(qname.c_str(), "SRR1171540.%u", &read);
    if (read > 1000000)
      continue;
    end_pos = 0;
    res[read] = CMAPPINGResult(rname, pos, end_pos, qname, mismatch, strand,
                               seq, qual);
  }
  fclose(fin);
}

void ReadBismarkToMRResults(const string& file_name,
                            const vector<CMAPPINGResult>& best_results,
                            const uint32_t& num_of_paried) {
  string line;
  string chrom;
  uint32_t start_pos;
  uint32_t end_pos;
  string read_name;
  int num_of_mismatches;
  char strand;
  string read_seq;
  string read_score;

  unsigned int read;
  char tmp[MAX_LINE_LENGTH];
  cerr << file_name << endl;
  ifstream fin(file_name);
  uint32_t cur_num_of_paried = 0;
  while (getline(fin, line)) {
    istringstream iss(line);
    iss >> chrom >> start_pos >> end_pos >> read_name >> num_of_mismatches
        >> strand >> read_seq >> read_score;
    if (read_name.substr(0, 5) == "FRAG:") {
      sscanf(read_name.c_str(), "FRAG:SRR1171540.%u_%s", &read, tmp);
      if (read > 1000000)
        continue;

      if (best_results[read].chrom == "XXX"
          && best_results[read].mismatch == 100) {
        continue;
      }

      if (num_of_mismatches == best_results[read].mismatch
          && chrom == best_results[read].chrom) {
        cur_num_of_paried++;
      } else {
        cout << "FRAG:SRR1171540." << read << endl;
        best_results[read].Output();
        cout << line << endl;
        cout << "---------------------------------------------------" << endl;
      }
    }
  }

  cout << cur_num_of_paried << "\t" << num_of_paried << "\t"
      << static_cast<double>(cur_num_of_paried) / num_of_paried << endl;

  fin.close();
}

int main(int argc, const char *argv[]) {
  vector<CMAPPINGResult> best_result(1000005);
  FILE * fin = fopen(argv[1], "r");
  char cline[MAX_LINE_LENGTH];
  char rname[MAX_LINE_LENGTH];
  char qname[MAX_LINE_LENGTH];
  char seq[MAX_LINE_LENGTH];
  char score[MAX_LINE_LENGTH];
  char strand;
  uint32_t start_pos, end_pos, mismatch;
  cerr << "read ground truth..." << endl;
  uint32_t read = 1;
  uint32_t num_of_paried = 0;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    sscanf(cline, "%s %u %u %s %u %c %s %s", rname, &start_pos, &end_pos, qname,
           &mismatch, &strand, seq, score);
    best_result[read] = CMAPPINGResult(rname, start_pos, end_pos, qname,
                                       mismatch, strand, seq, score);
    if (strcmp(rname, "XXX") != 0 && mismatch != 100) {
      num_of_paried++;
    }
    read++;
  }
  fclose(fin);

  cerr << "read mapping results..." << endl;
  vector<CMAPPINGResult> res(1000005);
  cerr << argv[3] << endl;

  if (strcmp(argv[2], "-bsmapper") == 0 || strcmp(argv[2], "-bsmap") == 0) {
    CompareMappingResultsMR(argv[3], best_result, num_of_paried);
  } else if (strcmp(argv[2], "-bismark") == 0) {
    ReadBismarkToMRResults(argv[3], best_result, num_of_paried);
  } else {
    cerr << "Please check the mapper..." << endl;
    return 0;
  }

  return 0;
}
