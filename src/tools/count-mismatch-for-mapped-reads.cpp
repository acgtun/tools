/*
 * count number of matched reads for different number of mismatches
 * unordered_map<uint32_t, uint32_t> mismatch_reads;
 * key: number of mismatches
 * value: number of matched reads
 */
#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "option.hpp"

#include "./../accuracy_bisulfite/read_sam.h"

using namespace std;

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);

  string mapping_file, mapper;
  bool paired = false;
  Option::GetOption("-f", mapping_file);
  Option::GetOption("-m", mapper);
  Option::ChkStrExist("-paired", paired);

  vector<CMAPPINGResult> res(50000005);
  cout << mapper << endl;
  cout << mapping_file << endl;
  cout << paired << endl;
  Read_SAM_Results(mapping_file.c_str(), res, mapper, paired);


  map<uint32_t, uint32_t> mismatch_count;
  for (uint32_t i = 0; i < res.size(); ++i) {

    uint32_t mismatch =
        paired ?  res[i].mismatch + res[i].mismatch2 : res[i].mismatch;
    //cout << res[i].mismatch <<  " " << res[i].mismatch2 << " " << res[i].mismatch + res[i].mismatch2 << " " << mismatch << endl;
    mismatch_count[mismatch]++;
  }
  cout << "size = " << res.size() << endl;
  for (map<uint32_t, uint32_t>::const_iterator it = mismatch_count.begin();
      it != mismatch_count.end(); ++it) {
    printf("%u\t%u\t%.2lf%%\n", it->first, it->second,
           100 * (double) it->second / (double) res.size());
  }

  return 0;
}
