#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "option.hpp"
#include "smithlab_os.hpp"

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

void CompareMappingResults(const string& file1, const string& file2,
                           ofstream& fout) {
  set<string> reads;
  map<string, CMAPPINGResult> res1, res2;
  ReadResult(file1, reads, res1);
  ReadResult(file2, reads, res2);

  uint32_t cnt_diff = 0, cnt_diff012 = 0, cnt_same012 = 0;
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
      if (ptr1->second.num_of_mismatches <= 2
          || ptr2->second.num_of_mismatches <= 2) {
        //printf("someting error~@!@@@\n");
        cnt_diff012++;
        /*
         * if (ptr1 != res1.end())
         cout << ptr1->second.read_name << endl;
         if (ptr2 != res2.end())
         cout << ptr2->second.read_name << endl;
         */
      }
      cnt_diff++;
      fout << "------------------------------------------------------" << endl;
    } else {
      if (ptr1->second.num_of_mismatches <= 2) {
        cnt_same012++;
      }
    }
  }


  cout << cnt_diff << "\t" << static_cast<double>(cnt_diff) / reads.size() << "\t";
  cout << cnt_diff012 << "\t" << static_cast<double>(cnt_diff012) / reads.size() << "\t";
  cout << cnt_same012  + cnt_diff012 << "\t" 
       << static_cast<double>(cnt_same012 + cnt_diff012) / reads.size() << "\t";
  cout  << static_cast<double>(cnt_diff) / reads.size() << "\t";
  cout << reads.size() << "\t";
  cout << "file2: " << file2 << endl;
}

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);

  string benchmark_file, mapping_file, output_file;
  Option::GetOption("-g", benchmark_file);
  Option::GetOption("-f", mapping_file);
  Option::GetOption("-o", output_file);

  vector<string> file_names;
  if (isdir(mapping_file.c_str())) {
    read_dir(mapping_file, "out", file_names);
  } else {
    file_names.push_back(mapping_file);
  }

  ofstream fout(output_file.c_str());
  cout << "diff\tprecentage\tcnt_diff012\tprecentage\treads has 012 mismatches\tprecentage\tmatched reads\tfile" << endl;
  for (uint32_t i = 0; i < file_names.size(); ++i) {
    CompareMappingResults(benchmark_file, file_names[i], fout);
  }
  fout.close();

  return 0;
}
