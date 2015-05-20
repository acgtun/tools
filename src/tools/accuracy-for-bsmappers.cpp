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
  CMAPPINGResult(string _chrom = "haoma", string _start_pos = "haoma",
                 string _end_pos = "haoma", string _read_name = "haoma",
                 int _num_of_mismatches = 0, char _strand = '+',
                 string _read_seq = "haoma", string _read_score = "haoma") {
    chrom = _chrom;
    start_pos = _start_pos;
    end_pos = _end_pos;
    read_name = _read_name;
    num_of_mismatches = _num_of_mismatches;
    strand = _strand;
    read_seq = _read_seq;
    read_score = _read_score;
  }
  friend bool CheckDifference(const CMAPPINGResult& r1,
                              const CMAPPINGResult& r2) {
    return r1.chrom == r2.chrom && r1.start_pos == r2.start_pos
        && r1.end_pos == r2.end_pos && r1.read_name == r2.read_name
        && r1.num_of_mismatches == r2.num_of_mismatches
        && r1.strand == r2.strand && r1.read_seq == r2.read_seq
        && r1.read_score == r2.read_score;
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

void ReadMRResult(const char* file_name, vector<CMAPPINGResult>& res) {
  string line;
  string chrom;
  string start_pos;
  string end_pos;
  string read_name;
  int num_of_mismatches;
  char strand;
  string read_seq;
  string read_score;

  unsigned int read;
  cerr << file_name << endl;
  ifstream fin(file_name);
  while (getline(fin, line)) {
    istringstream iss(line);
    iss >> chrom >> start_pos >> end_pos >> read_name >> num_of_mismatches
        >> strand >> read_seq >> read_score;
    sscanf(read_name.c_str(), "SRR1171540.%u", &read);
    if (read > 1000000)
      continue;
    res[read] = CMAPPINGResult(chrom, start_pos, end_pos, read_name,
                               num_of_mismatches, strand, read_seq, read_score);
  }
  fin.close();
}

void ReadBsmapResults(const char* file_name, vector<CMAPPINGResult>& res) {
  FILE * fin = fopen(file_name, "r");
  if (!fin) {
    cerr << "cannot open input file " << file_name << endl;
  }
  string qname, rname, cigar, rnext, seq, qual, pos, end_pos, NM, MD;
  uint32_t flag, mapq, pnext, tlen;
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
    res[read] = CMAPPINGResult(rname, pos, end_pos, qname, mismatch, strand,
                               seq, qual);
  }
  fclose(fin);
}

void ReadBismarkResults(const char* file_name, vector<CMAPPINGResult>& res) {
  FILE * fin = fopen(file_name, "r");
  if (!fin) {
    cerr << "cannot open input file " << file_name << endl;
  }
  string qname, rname, cigar, rnext, seq, qual, pos, end_pos, NM, MD;
  uint32_t flag, mapq, pnext, tlen;
  int mismatch;
  char strand = '+';
  char cline[MAX_LINE_LENGTH];
  char tmp[MAX_LINE_LENGTH];
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    if (cline[0] == '@')
      continue;
    istringstream iss(cline);
    cerr << cline << endl;
    iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext
        >> tlen >> seq >> qual >> NM >> MD;
    sscanf(NM.c_str(), "NM:i:%d", &mismatch);
    unsigned int read;
    sscanf(qname.c_str(), "SRR1171540.%u_%s", &read, tmp);
    if (read > 1000000)
      continue;
    cerr << read << " " << rname << " " << mismatch << endl;
    res[read] = CMAPPINGResult(rname, pos, end_pos, qname, mismatch, strand,
                               seq, qual);
  }
  fclose(fin);
}

void ReadBismarkToMRResults(const char* file_name,
                            vector<CMAPPINGResult>& res) {
  string line;
  string chrom;
  string start_pos;
  string end_pos;
  string read_name;
  int num_of_mismatches;
  char strand;
  string read_seq;
  string read_score;

  unsigned int read;
  char tmp[MAX_LINE_LENGTH];
  cerr << file_name << endl;
  ifstream fin(file_name);
  while (getline(fin, line)) {
    istringstream iss(line);
    iss >> chrom >> start_pos >> end_pos >> read_name >> num_of_mismatches
        >> strand >> read_seq >> read_score;
    sscanf(read_name.c_str(), "SRR1171540.%u_%s", &read, tmp);
    if (read > 1000000)
      continue;
    res[read] = CMAPPINGResult(chrom, start_pos, end_pos, read_name,
                               num_of_mismatches, strand, read_seq, read_score);
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
  ofstream fout("err.see.txt");
  uint32_t correct = 0;
  for (uint32_t i = 1; i <= 1000000; ++i) {
    int uni = unique(count[i]);
    if (uni == -1 || uni == MAX_INT) {
      if (res[i].chrom == "haoma") {
        correct++;
      } else {
        fout << "duo le: " << i << endl;
        res[i].Output(fout);
      }
    } else {
      if (uni == res[i].num_of_mismatches) {
        correct++;
      } else {
        if (res[i].chrom == "haoma") {
          fout << "mei zhao dao" << endl;
        } else {
          fout << "number of mismatch not identical: " << i << endl;
          fout << count[i][0] << " " << count[i][1] << " " << count[i][2] << " "
               << count[i][3] << " " << count[i][4] << " " << count[i][5] << " "
               << count[i][6] << endl;
          res[i].Output(fout);
        }
      }
    }
  }

  cout << correct << "\t" << static_cast<double>(correct) / 1000000 << endl;
}

int main(int argc, const char *argv[]) {
  vector<vector<uint32_t> > count(1000005, vector<uint32_t>(7, 0));
  FILE * fin = fopen(argv[1], "r");
  char cline[MAX_LINE_LENGTH];
  cerr << "read ground truth..." << endl;
  unsigned int read, m0, m1, m2, m3, m4, m5, m6;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    sscanf(cline, "SRR1171540.%u %u %u %u %u %u %u %u", &read, &m0, &m1, &m2,
           &m3, &m4, &m5, &m6);
    count[read][0] = m0;
    count[read][1] = m1;
    count[read][2] = m2;
    count[read][3] = m3;
    count[read][4] = m4;
    count[read][5] = m5;
    count[read][6] = m6;
  }
  fclose(fin);

  return 0;

  cerr << "read mapping results..." << endl;
  vector<CMAPPINGResult> res(1000005);
  cerr << argv[3] << endl;
  if (strcmp(argv[2], "-bsmap") == 0) {
    ReadBsmapResults(argv[3], res);
  } else if (strcmp(argv[2], "-bismark") == 0) {
    ReadBismarkToMRResults(argv[3], res);
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
