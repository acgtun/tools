#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <cstdint>

using namespace std;

int main(int argc, const char *argv[]) {
  /* input file is methylation level for each CpG site (symmetric) */

  string chrom, CpG;
  uint32_t position;
  char strand;
  double meth_lev;
  uint32_t coverage;

  string outfile = argv[1];
  outfile += "_meth_lev";
  ifstream fin(argv[1]);
  ofstream fout(outfile.c_str());
  cout << "input " << argv[1] << endl;
  cout << "output " << outfile << endl;

  map<string, map<uint32_t, pair<double, uint32_t> > > map_meth_lve;

  uint32_t total_CpG_pos = 0;
  uint32_t coverage_less_2 = 0;
  uint32_t coverage_less_5 = 0;
  uint32_t coverage_less_8 = 0;
  uint32_t coverage_less_10 = 0;
  while (fin >> chrom >> position >> strand >> CpG >> meth_lev >> coverage) {
    total_CpG_pos++;
    if(coverage < 2) coverage_less_2++;
    if(coverage < 5) coverage_less_5++;
    if(coverage < 8) coverage_less_8++;
    if(coverage < 10) coverage_less_10++;

    map_meth_lve[chrom][position] = make_pair(meth_lev, coverage);
  }
  cout << coverage_less_2 << " " << total_CpG_pos << " " << coverage_less_2 / (double)total_CpG_pos << endl;
  cout << coverage_less_5 << " " << total_CpG_pos << " " << coverage_less_5 / (double)total_CpG_pos << endl;
  cout << coverage_less_8 << " " << total_CpG_pos << " " << coverage_less_8 / (double)total_CpG_pos << endl;
  cout << coverage_less_10 << " " << total_CpG_pos << " " << coverage_less_10 / (double)total_CpG_pos << endl;


  system("mkdir -p ./chrom_meth_lve/");
  for (map<string, map<uint32_t, pair<double, uint32_t> > >::iterator it = map_meth_lve.begin();
      it != map_meth_lve.end(); ++it) {
    string chrom_outfile = "./chrom_meth_lve/";
    chrom_outfile += it->first + "_chrom_meth_lev";
    ofstream fchrom(chrom_outfile.c_str());
    for (map<uint32_t, pair<double, uint32_t> >::iterator it2 = it->second.begin();
        it2 != it->second.end(); ++it2) {
      //fout << it->first << "\t" << it2->first << "\t" << it2->second << endl;
      fout << it2->second.first << "\t" << it2->second.second << endl;
      fchrom << it2->first << "\t" << it2->second.first << endl;
    }
    fchrom.close();
  }

  fin.close();
  fout.close();

  return 0;
}
