#include <map>
#include <set>
#include <vector>
#include <string>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

struct MethSite {
  MethSite(const string& _chrom, const uint32_t& _position,
           const double& _meth_lev, const int& _coverage)
      : chrom(_chrom),
        position(_position),
        meth_lev(_meth_lev),
        coverage(_coverage) {
  }
  string chrom;
  uint32_t position;
  double meth_lev;
  int coverage;
};

void ReadMethCountsResults(const string& file, vector<MethSite>& lev) {
  cout << "file: " << file << endl;
  string chrom, CpG;
  uint32_t position;
  char strand;
  double meth_lev;
  int coverage;

  ifstream fin(file.c_str());

  while (fin >> chrom >> position >> strand >> CpG >> meth_lev >> coverage) {
    lev.push_back(MethSite(chrom, position, meth_lev, coverage));
  }

  fin.close();
}

bool check(const MethSite& site1, const MethSite& site2) {
  //cout << site1.chrom << "\t" << site1.position << endl;
  //cout << site2.chrom << "\t" << site2.position << endl;
  return site1.chrom == site2.chrom && site1.position == site2.position;
}

int main(int argc, const char *argv[]) {
  /* input file is methylation level for each CpG site (symmetric) */
  vector<MethSite> lev_1;
  vector<MethSite> lev_2;

  ReadMethCountsResults(argv[1], lev_1);
  ReadMethCountsResults(argv[2], lev_2);

  string file1 = argv[1];
  string file2 = argv[2];

  size_t p1 = file1.find_last_of('/');
  size_t p2 = file2.find_last_of('/');
  file1[p1] = '_';
  file2[p2] = '_';
  //file1 = file1.substr(p1 + 1);
  //file2 = file2.substr(p2 + 1);
  cout << p1 << " " << p2 << endl;

  if (lev_1.size() != lev_2.size()) {
    fprintf(stderr, "the size should be the same\n");
    exit(0);
  }

  for (int d = 2; d <= 10; ++d) {
    char file_name[100];
    sprintf(file_name, "%s_vs_%s_diff%d_%s.meth_lev_comp", file1.c_str(),
            file2.c_str(), d, lev_1[0].chrom.c_str());
    ofstream fout(file_name);
    cout << "output " << file_name << endl;

    for (size_t i = 0; i < lev_1.size(); ++i) {
      if (i != 0 && lev_1[i].chrom != lev_2[i - 1].chrom) {
        fout.close();
        char file_name[100];
        sprintf(file_name, "%s_vs_%s_diff%d_%s.meth_lev_comp", file1.c_str(),
                file2.c_str(), d, lev_1[i].chrom.c_str());
        fout.open(file_name);
        cout << "output " << file_name << endl;
      }
      if (!check(lev_1[i], lev_2[i])) {
        fprintf(stderr, "positions should be the same\n");
        cout << lev_1[i].chrom << "\t" << lev_1[i].position << endl;
        cout << lev_2[i].chrom << "\t" << lev_2[i].position << endl;
        exit(0);
      } else {
        if (abs(lev_1[i].coverage - lev_2[i].coverage) >= d && fabs(lev_1[i].meth_lev - lev_2[i].meth_lev) > 0.4)
//          fout << lev_1[i].chrom << "\t" << lev_1[i].position << "\t"
//               << lev_1[i].coverage << "\t" << lev_2[i].coverage << endl;
          fout << i << "\t" << lev_1[i].position << "\t" << lev_1[i].coverage
               << "\t" << lev_2[i].coverage << endl;
      }
    }

    fout.close();
  }

  return 0;
}
