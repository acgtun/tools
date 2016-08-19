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

#include "read_sam.h"

using namespace std;

int main(int argc, const char *argv[]) {
  cerr << "read mapping results..." << endl;
  InitProgram(argc, argv);
  int num_of_reads = 50000005;
  bool paired = false;
  Option::GetOption("-s", num_of_reads, 50000005);
  Option::ChkStrExist("-paired", paired);

  vector<CMAPPINGResult> res(num_of_reads);
  Read_SAM_Results(argv[2], res, argv[1], paired);

  map<int, int> count_mismatch;
  for (size_t i = 0; i < res.size(); ++i) {
    int num_of_mismatches = 0;
    num_of_mismatches = res[i].mismatch;
    if (paired) {
      num_of_mismatches += res[i].mismatch2;
    }
    count_mismatch[num_of_mismatches]++;
  }

  for (map<int, int>::iterator it = count_mismatch.begin();
      it != count_mismatch.end(); ++it) {
    cout << it->first << " " << it->second << endl;
  }

  return 0;
}
