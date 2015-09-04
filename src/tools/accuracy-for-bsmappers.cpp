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

using namespace std;

const int MAX_LINE_LENGTH = 10000;
const int MAX_INT = std::numeric_limits<int>::max();

struct CMAPPINGResult {
  CMAPPINGResult(string _chrom = "-1", string _start_pos = "-1",
                 string _end_pos = "-1", string _read_name = "-1",
                 int _mismatches = -1, char _strand = '+',
                 string _read_seq = "-1", string _read_score = "-1") {
    chrom = _chrom;
    start_pos = _start_pos;
    end_pos = _end_pos;
    read_name = _read_name;
    mismatches = _mismatches;
    strand = _strand;
    read_seq = _read_seq;
    read_score = _read_score;
  }
  friend bool CheckDifference(const CMAPPINGResult& r1,
                              const CMAPPINGResult& r2) {
    return r1.chrom == r2.chrom && r1.start_pos == r2.start_pos
        && r1.end_pos == r2.end_pos && r1.read_name == r2.read_name
        && r1.mismatches == r2.mismatches
        && r1.strand == r2.strand && r1.read_seq == r2.read_seq
        && r1.read_score == r2.read_score;
  }

  void Output(ofstream& fout) {
    fout << chrom << "\t" << start_pos << "\t" << end_pos << "\t" << read_name
         << "\t" << mismatches << "\t" << strand << "\t" << read_seq
         << "\t" << read_score << endl;
  }
  string chrom;
  string start_pos;
  string end_pos;
  string read_name;
  int mismatches;
  char strand;
  string read_seq;
  string read_score;
};

void ReadMRResult(const char* file_name, vector<CMAPPINGResult>& res) {
  string line;
  string chrom;
  string start_pos;
  string end_pos;
  string read_name;
  int mismatches;
  char strand;
  string read_seq;
  string read_score;

  unsigned int read;
  cerr << file_name << endl;
  unsigned int SRRName, line_count = 0;
  ifstream fin(file_name);
  while (getline(fin, line)) {
    istringstream iss(line);
    iss >> chrom >> start_pos >> end_pos >> read_name >> mismatches
        >> strand >> read_seq >> read_score;
    cout << read_name << endl;
    sscanf(read_name.c_str(), "SRR%u.%u", &SRRName, &read);
    if(line_count < 10) {
      cout << SRRName << endl;
      cout << line << endl;
    }
    line_count++;
    if (read > 1000000)
      break;
    res[read] = CMAPPINGResult(chrom, start_pos, end_pos, read_name,
                               mismatches, strand, read_seq, read_score);
  }
  fin.close();
}


int unique(const vector<uint32_t>& mismatch) {
  for (int i = 0; i <= 6; ++i) {
    if (mismatch[i] == 0)
      continue;
    if (mismatch[i] == 1)
      return i;
    if (mismatch[i] >= 2)
      return MAX_INT;
  }
  return -1;
}

void CompareMappingResults(vector<CMAPPINGResult>& res,
                           const vector<vector<uint32_t> >& count) {
  cout << "comparing...." << endl;
  int TP = 0, FP = 0, FN = 0, TN = 0;
  int tp[7] = { 0 }, fp[7] = { 0 }, fn[7] = { 0 };
  for (uint32_t i = 1; i <= 1000000; ++i) {
    int uni = unique(count[i]);
    if (uni == -1 || uni == MAX_INT) {
      if (res[i].chrom == "-1") {
        TN++;
        //tn[res[i].mismatches]++;
      } else {
        FP++;
        if (res[i].mismatches > 6) {
          fp[6]++;
        } else {
          fp[res[i].mismatches]++;
        }
      }
    } else {
      if (uni == res[i].mismatches) {
        TP++;
        tp[uni]++;
      } else {
        FN++;
        fn[uni]++;
      }
    }
  }

  printf("TP: %d TN:%d FP:%d FN:%d\n", TP, TN, FP, FN);
  printf("Total:     %.4lf %.4lf %.4lf\n", (double) TP / (double) (TP + FN),
         (double) TP / (double) (TP + FP),
         (double) 2.0 * TP / (double) (2.0 * TP + FP + FN));

  for (int i = 0; i <= 6; ++i) {
    printf("%.4lf %.4lf %.4lf\n",
           (double) tp[i] / (double) (tp[i] + fn[i]),
           (double) tp[i] / (double) (tp[i] + fp[i]),
           (double) 2.0 * tp[i] / (double) (2.0 * tp[i] + fp[i] + fn[i]));
  }
}

int main(int argc, const char *argv[]) {
  /* input summary of count positions on diff mismatch*/
  vector<vector<uint32_t> > count(1000005, vector<uint32_t>(7, 0));
  FILE * fin = fopen(argv[1], "r");
  char cline[MAX_LINE_LENGTH];
  cerr << "read ground truth..." << endl;
  unsigned int SRRName;
  unsigned int line_count = 0;
  unsigned int read, m0, m1, m2, m3, m4, m5, m6;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    sscanf(cline, "SRR%u.%u %u %u %u %u %u %u %u", &SRRName, &read, &m0, &m1, &m2,
           &m3, &m4, &m5, &m6);
    if(line_count < 10) {
      cout << SRRName << endl;
      cout << cline << endl;
    }
    line_count++;
    count[read][0] = m0;
    count[read][1] = m1;
    count[read][2] = m2;
    count[read][3] = m3;
    count[read][4] = m4;
    count[read][5] = m5;
    count[read][6] = m6;
  }
  fclose(fin);

  cerr << "read mapping results..." << endl;
  vector<CMAPPINGResult> res(1000005);
  cerr << argv[3] << endl;
  if (strcmp(argv[2], "-bsmap") == 0) {
    ReadMRResult(argv[3], res);
  } else if (strcmp(argv[2], "-bismark") == 0) {
    ReadMRResult(argv[3], res);
  } else if (strcmp(argv[2], "-bsmapper") == 0) {
    ReadMRResult(argv[3], res);
  } else {
    cerr << "Please check the mapper..." << endl;
    return 0;
  }

  cerr << "comparing..." << endl;
  CompareMappingResults(res, count);

  return 0;
}
