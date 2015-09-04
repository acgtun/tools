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
  CMAPPINGResult(string _chrom = "-1", uint32_t _start_pos = 0,
                 uint32_t _end_pos = 0, string _read_name = "-1",
                 int _mismatch = -1, char _strand = '+',
                 string _read_seq = "-1", string _read_score = "-1") {
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

//SRR1287225
//SRR1171540
void CompareMappingResultsMR(const string& file_name,
                             const vector<CMAPPINGResult>& best_results,
                             const uint32_t& num_of_paried,
                             map<int, int>& count_on_mismatch, const int& sum_paired) {
  cerr << "CompareMappingResultsMR..." << endl;
  //cerr << file_name << endl;
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
  cerr << file_name << endl;
  int TP = 0, FP = 0, FN = 0, TN = 0;
  int tp[7] = { 0 }, fp[7] = { 0 }, fn[7] = { 0 };
  ifstream fin(file_name);
  set<unsigned int> read_set;
  unsigned int SRRName;
  while (getline(fin, line)) {
    istringstream iss(line);
    iss >> chrom >> start_pos >> end_pos >> read_name >> num_of_mismatches
        >> strand >> read_seq >> read_score;
    //cout << "start" << endl;
    //cout << line << endl;
    bool paired = false;
    if (read_name.substr(0, 5) == "FRAG:") {
      paired = true;
      sscanf(read_name.c_str(), "FRAG:SRR%u.%u", &SRRName, &read);
    } else {
      sscanf(read_name.c_str(), "SRR%u.%u", &SRRName, &read);
    }
    if (read > 1000000)
      break;
    if (read_set.find(read) != read_set.end()) {
      continue;
    }

    read_set.insert(read);

    if (best_results[read].chrom == "XXX"
        && best_results[read].mismatch == 100) {
      if (paired) {
        FP++;
        if (num_of_mismatches > 6) {
          fp[6]++;
        } else {
          fp[num_of_mismatches]++;
        }
      } else {
        TN++;
      }
    } else {
      if (paired && num_of_mismatches == best_results[read].mismatch
          && chrom == best_results[read].chrom) {
        TP++;
        tp[num_of_mismatches]++;
      } else {
        FN++;
        if (num_of_mismatches > 6) {
          fn[6]++;
        } else {
          fn[num_of_mismatches]++;
        }
      }
    }
  }

  printf("TP: %d TN:%d FP:%d FN:%d\n", TP, TN, FP, sum_paired - TP);
  cout << "TP + FN = " << TP + FN << endl;
  FN = sum_paired - TP;
  printf("Total:     %.4lf %.4lf %.4lf\n", (double) TP / (TP + FN),
         (double) TP / (double) (TP + FP),
         (double) 2.0 * TP / (double) (2.0 * TP + FP + FN));

  for (int i = 0; i <= 6; ++i) {
    printf("%d: tp: %d fp:%d fn:%d\n", i, tp[i], fp[i], count_on_mismatch[i] - tp[i]);
  }
  printf("----------------------\n");
  for (int i = 0; i <= 6; ++i) {
    //printf("tp: %d fp:%d fn:%d\n", tp[i],  fp[i], fn[i]);
    fn[i] = count_on_mismatch[i] - tp[i];
    printf("%.4lf %.4lf %.4lf\n",
           (double) tp[i] / (double) (tp[i] + fn[i]),
           (double) tp[i] / (double) (tp[i] + fp[i]),
           (double) 2.0 * tp[i] / (double) (2.0 * tp[i] + fp[i] + fn[i]));
  }

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
  unsigned int SRRName, line_count = 0;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    if (cline[0] == '@')
      continue;
    istringstream iss(cline);
    iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext
        >> tlen >> seq >> qual >> NM >> MD;
    sscanf(NM.c_str(), "NM:i:%d", &mismatch);
    unsigned int read;
    sscanf(qname.c_str(), "SRR%u.%u", &SRRName, &read);
    if (line_count < 10) {
      //cout << SRRName << endl;
      //cout << cline << endl;
    }
    line_count++;
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
  unsigned int SRRName, line_count = 0;
  cerr << file_name << endl;
  ifstream fin(file_name);
  uint32_t cur_num_of_paried = 0;
  while (getline(fin, line)) {
    istringstream iss(line);
    iss >> chrom >> start_pos >> end_pos >> read_name >> num_of_mismatches
        >> strand >> read_seq >> read_score;
    if (read_name.substr(0, 5) == "FRAG:") {
      sscanf(read_name.c_str(), "FRAG:SRR%u.%u_%s", &SRRName, &read, tmp);
      if (line_count < 10) {
        //cout << SRRName << endl;
        //cout << line << endl;
      }
      line_count++;
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
        //cout << "FRAG:" << SRRName << "." << read << endl;
        //best_results[read].Output();
        //cout << line << endl;
        //cout << "---------------------------------------------------" << endl;
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
  map<int, int> count_on_mismatch;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    sscanf(cline, "%s %u %u %s %u %c %s %s", rname, &start_pos, &end_pos, qname,
           &mismatch, &strand, seq, score);
    best_result[read] = CMAPPINGResult(rname, start_pos, end_pos, qname,
                                       mismatch, strand, seq, score);

    if (strcmp(rname, "XXX") != 0 && mismatch != 100) {
      num_of_paried++;
      count_on_mismatch[mismatch]++;
    }
    read++;
  }
  int sum = 0;
  for (map<int, int>::iterator it = count_on_mismatch.begin();
      it != count_on_mismatch.end(); ++it) {
    cout << it->first << " " << it->second << endl;
    sum += it->second;
  }
  cout << "sum = " << sum << endl;
  fclose(fin);

  cerr << "read mapping results..." << endl;
  vector<CMAPPINGResult> res(1000005);
  cerr << argv[3] << endl;

  if (strcmp(argv[2], "-bsmapper") == 0 || strcmp(argv[2], "-bsmap") == 0) {
    CompareMappingResultsMR(argv[3], best_result, num_of_paried, count_on_mismatch, sum);
  } else if (strcmp(argv[2], "-bismark") == 0) {
    ReadBismarkToMRResults(argv[3], best_result, num_of_paried);
  } else {
    cerr << "Please check the mapper..." << endl;
    return 0;
  }

  return 0;
}
