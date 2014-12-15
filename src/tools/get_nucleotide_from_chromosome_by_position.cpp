#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

int main(int argc, const char *argv[]) {
  ifstream fin(argv[1]);
  string line;
  string chrom_seq;
  while (getline(fin, line)) {
    if (line[0] == '>')
      continue;
    else {
      chrom_seq += line;
    }
  }
  uint64_t start;
  istringstream iss(argv[2]);
  iss >> start;
  cout << "size of chromosome " << argv[1] << " is " << chrom_seq.size() << endl;
  cout << "start position: " << start << endl;
  for (uint64_t i = start, j = 0; i < chrom_seq.size() && j < 100; ++i, ++j) {
    cout << chrom_seq[i];
  }
  cout << endl;
  return 0;
}
