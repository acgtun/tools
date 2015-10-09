#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

const int MAX_LINE_LENGTH = 1000;

struct MatchResult {
  MatchResult(string _chrom = "XX", uint32_t _pos = 0, char _strand = '#',
              uint32_t _mismatch = 100)
      : chrom(_chrom),
        pos(_pos),
        strand(_strand),
        mismatch(_mismatch) {
  }

  string chrom;
  uint32_t pos;
  char strand;
  uint32_t mismatch;
};

void GetCounts(const string& count_file, vector<vector<uint32_t> >& count,
               vector<uint32_t>& threshold) {
  /* input summary of count positions on diff mismatch*/
  cout << count_file << endl;
  FILE * fin = fopen(count_file.c_str(), "r");
  char cline[MAX_LINE_LENGTH];
  cerr << "read ground truth..." << endl;
  unsigned int SRRName, line_count = 0;
  unsigned int read, m0, m1, m2, m3, m4, m5, m6;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    //cerr << cline << endl;
    //cout << "---------------------------------------" << endl;
    //cout << cline << endl;
    sscanf(cline, "SRR%u.%u%u%u%u%u%u%u%u", &SRRName, &read, &m0, &m1, &m2, &m3, &m4,
           &m5, &m6);
    //if(line_count < 10) {
      //cout << SRRName << endl;
      //cout << "read:" << read << endl;
      //cout << cline << endl;
    //}
    line_count++;
    //cerr << read << " " << m0 << " " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << m5 << " " << m6 << endl;
    //cout << "---------------------------------------" << endl;
    count[read][0] = m0;
    count[read][1] = m1;
    count[read][2] = m2;
    count[read][3] = m3;
    count[read][4] = m4;
    count[read][5] = m5;
    count[read][6] = m6;
  }
  fclose(fin);
  cout << "hehhh" << endl;
  for (uint32_t i = 1; i <= 1000000; ++i) {
    uint32_t sum = 0;
    for (uint32_t j = 0; j <= 6; ++j) {
      sum += count[i][j];
      if (sum >= 1000) {
        threshold[i] = j;
        break;
      }
    }
    if (sum < 1000) {
      threshold[i] = 6;
    }
  }
  string threshold_output = count_file + "_threshold.txt";
  ofstream fout(threshold_output.c_str());
  for (uint32_t i = 1; i <= 1000000; ++i) {
    fout << i << " " << threshold[i] << endl;
  }
  fout.close();
}

int main(int argc, const char **argv) {
  // input 1  .count.txt.sum.txt
  // input 2 .._post_processing.txt.count.txt
  vector<uint32_t> threshold(1000005, 0);
  vector<vector<uint32_t> > count(1000005, vector<uint32_t>(7, 0));
  vector<vector<MatchResult> > match_results(1000005,
                                             vector<MatchResult>(1000));
  vector<uint32_t> match_results_size(1000005, 0);

  GetCounts(argv[1], count, threshold);
  FILE * fin = fopen(argv[2], "r");
  char cline[MAX_LINE_LENGTH];
  char chrom[MAX_LINE_LENGTH];
  uint32_t pos;
  char strand;
  uint32_t mismatch;
  uint32_t read;
  unsigned int SRRName, line_count = 0;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    sscanf(cline, "SRR%u.%u %s %u %c %u", &SRRName, &read, chrom, &pos, &strand,
           &mismatch);
    if(line_count < 10) {
      cout << SRRName << endl;
      cout << cline << endl;
    }
    line_count++;
    if (mismatch <= threshold[read] && match_results_size[read] < 1000) {
      match_results[read][match_results_size[read]] = MatchResult(chrom, pos,
                                                                  strand,
                                                                  mismatch);
      match_results_size[read]++;
    }
  }
  fclose(fin);

  string output = argv[2];
  output += "top_1000.txt";
  FILE * fout = fopen(output.c_str(), "w");
  for (uint32_t i = 1; i <= 1000000; ++i) {
    fprintf(fout, "%u %u\n", i, match_results_size[i]);
    for (uint32_t j = 0; j < match_results_size[i]; ++j) {
      fprintf(fout, "%s %u %c %u\n", match_results[i][j].chrom.c_str(),
              match_results[i][j].pos, match_results[i][j].strand,
              match_results[i][j].mismatch);
    }
  }
  fclose(fout);

  return 0;
}
