#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

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

void ReadResult(const string& file_name, set<string>& reads,
                map<string, CMAPPINGResult>& res) {
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
    reads.insert(read_name);
    res.insert(
        make_pair(
            read_name,
            CMAPPINGResult(chrom, start_pos, end_pos, read_name,
                           num_of_mismatches, strand, read_seq, read_score)));
  }
}

int main(int argc, const char *argv[]) {

  map<string, CMAPPINGResult> res1, res2;
  set<string> reads;

  ReadResult(argv[1], reads, res1);
  ReadResult(argv[2], reads, res2);

  ofstream fout("diff.txt");
  for (set<string>::const_iterator it = reads.begin(); it != reads.end();
      ++it) {
    map<string, CMAPPINGResult>::iterator ptr1 = res1.find(*it);
    map<string, CMAPPINGResult>::iterator ptr2 = res2.find(*it);
    if (ptr1 == res1.end() || ptr2 == res2.end()
        || !CheckDifference(ptr1->second, ptr2->second)) {
      if (ptr1 != res1.end()) {
        fout << "1@ ";
        ptr1->second.Output(fout);
      }
      if (ptr2 != res2.end()) {
        fout << "2@ ";
        ptr2->second.Output(fout);
      }
      fout << "------------------------------------------------------" << endl;
    }
  }
  fout.close();

  return 0;
}
