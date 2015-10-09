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

#include "read_sam.h"

void CompareMappingResults(vector<CMAPPINGResult>& res,
                           vector<CMAPPINGResult>& best_result) {
  cout << "comparing...." << endl;
  int TP = 0, FP = 0, FN = 0, TN = 0;
  int tp[7] = { 0 }, fp[7] = { 0 }, fn[7] = { 0 };
  for (uint32_t i = 1; i <= 1000000; ++i) {
    if (best_result[i].chrom == "XXX" && best_result[i].mismatch == 100) {
      if (res[i].chrom != "XXX") {
        FP++;
        if (res[i].mismatch > 6) {
          fp[6]++;
        } else {
          fp[res[i].mismatch]++;
        }
      } else {
        TN++;
      }
    } else {
      if (best_result[i].chrom == res[i].chrom
          && best_result[i].start_pos + 1 == res[i].start_pos) {
        TP++;
        tp[best_result[i].mismatch]++;
      } else {
        FN++;
        fn[best_result[i].mismatch]++;
      }
    }
  }

  fprintf(stderr, "TP: %d TN:%d FP:%d FN:%d\n", TP, TN, FP, FN);
  fprintf(stderr, "Total:\n%.4lf\n%.4lf\n%.4lf\n",
          (double) TP / (double) (TP + FN), (double) TP / (double) (TP + FP),
          (double) 2.0 * TP / (double) (2.0 * TP + FP + FN));
  fprintf(stderr, "-------------\n");
  for (int i = 0; i <= 6; ++i) {
    fprintf(stderr, "%.4lf\n",
            (double) 2.0 * tp[i] / (double) (2.0 * tp[i] + fp[i] + fn[i]));
  }
}

int main(int argc, const char *argv[]) {
  /* input summary of count positions on diff mismatch*/
  vector<CMAPPINGResult> best_result(1000005);
  FILE * fin = fopen(argv[1], "r");
  char cline[MAX_LINE_LENGTH];
  char rname[MAX_LINE_LENGTH];
  unsigned int line_count = 0;
  uint32_t read = 1;
  uint32_t start_pos, mismatch;
  map<int, int> count_on_mismatch;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    sscanf(cline, "%s %u %u", rname, &start_pos, &mismatch);
    if (line_count < 1) {
      //cout << SRRName << endl;
      cout << cline << endl;
    }
    line_count++;
    best_result[read] = CMAPPINGResult(rname, start_pos, mismatch, 0, 100);
    if (strcmp(rname, "XXX") != 0 && mismatch != 100) {
      count_on_mismatch[mismatch]++;
    }
    read++;
  }
  fclose(fin);

  int sum = 0;
  for (map<int, int>::iterator it = count_on_mismatch.begin();
      it != count_on_mismatch.end(); ++it) {
    cout << it->first << " " << it->second << endl;
    sum += it->second;
  }
  cout << "sum = " << sum << endl;

  cerr << "read mapping results..." << endl;
  vector<CMAPPINGResult> res(1000005);
  cerr << argv[3] << endl;
  Read_SAM_Results(argv[3], res, argv[2]);

  CompareMappingResults(res, best_result);

  return 0;
}
