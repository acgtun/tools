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
#include "read_brat_bw.h"

using namespace std;

void CompareMappingResultsMR(const vector<CMAPPINGResult>& res,
                             const vector<CMAPPINGResult>& best_results,
                             const uint32_t& num_of_paried,
                             map<int, int>& count_on_mismatch,
                             const int& sum_paired) {
  cerr << "CompareMappingResultsMR..." << endl;
  int num_of_mismatches;
  int TP = 0, FP = 0, FN = 0, TN = 0;
  int tp[7] = { 0 }, fp[7] = { 0 }, fn[7] = { 0 };
  for (int i = 1; i <= 1000000; ++i) {
    cout << res[i].start_pos << " " << res[i].start_pos2 << " , "
        << best_results[i].start_pos << " " << best_results[i].start_pos2
        << endl;

    int paired = 0;
    if (res[i].start_pos != 0 && res[i].start_pos2 != 0) {
      paired = 1;
    }

    /////////////////////////////////////////////
    //////////////////////////////////////////
    if (best_results[i].chrom == "XXX") {
      if (paired) {
        FP++;
        num_of_mismatches = res[i].mismatch + res[i].mismatch2;
        if (num_of_mismatches > 6) {
          fp[6]++;
        } else {
          fp[num_of_mismatches]++;
        }
      } else {
        TN++;
      }
    } else {
      num_of_mismatches = best_results[i].mismatch + best_results[i].mismatch2;
      /*if ((res[i].start_pos > best_results[i].start_pos
       && res[i].start_pos - best_results[i].start_pos < 5)
       || (res[i].start_pos < best_results[i].start_pos
       && best_results[i].start_pos - res[i].start_pos < 100005)) {
       cout << "best " << res[i].start_pos << " " << best_results[i].start_pos
       << " " << res[i].start_pos2 << " " << best_results[i].start_pos2
       << endl;
       }*/
      /////////////////////////////////
      if (paired && res[i].start_pos == best_results[i].start_pos + 1
          /*&& res[i].start_pos2 == best_results[i].start_pos2 + 1*/
          && res[i].chrom == best_results[i].chrom) {
        TP++;
        tp[num_of_mismatches]++;
      } else {
        FN++;
        if (num_of_mismatches > 6) {
          fn[6]++;
        } else {
          fn[num_of_mismatches]++;
        }
      }
    }
    printf("TP: %d TN:%d FP:%d FN:%d\n", TP, TN, FP, sum_paired - TP);
  }

  printf("TP: %d TN:%d FP:%d FN:%d\n", TP, TN, FP, sum_paired - TP);
  cout << "TP + FN = " << TP + FN << endl;
  FN = sum_paired - TP;
  printf("Total:\n%.4lf\n%.4lf\n%.4lf\n", (double) TP / (TP + FN),
         (double) TP / (double) (TP + FP),
         (double) 2.0 * TP / (double) (2.0 * TP + FP + FN));
  fprintf(stderr, "-------------\n");
  for (int i = 0; i <= 6; ++i) {
    printf("%d: tp: %d fp:%d fn:%d\n", i, tp[i], fp[i],
           count_on_mismatch[i] - tp[i]);
  }
  printf("----------------------\n");
 /*for (int i = 0; i <= 6; ++i) {
    //printf("tp: %d fp:%d fn:%d\n", tp[i],  fp[i], fn[i]);
    fn[i] = count_on_mismatch[i] - tp[i];
    printf("%.4lf %.4lf %.4lf\n", (double) tp[i] / (double) (tp[i] + fn[i]),
           (double) tp[i] / (double) (tp[i] + fp[i]),
           (double) 2.0 * tp[i] / (double) (2.0 * tp[i] + fp[i] + fn[i]));
  }*/

  for (int i = 0; i <= 6; ++i) {
    //printf("tp: %d fp:%d fn:%d\n", tp[i],  fp[i], fn[i]);
    fn[i] = count_on_mismatch[i] - tp[i];
    printf("%.4lf\n",
           (double) 2.0 * tp[i] / (double) (2.0 * tp[i] + fp[i] + fn[i]));
  }


}

int main(int argc, const char *argv[]) {
  vector<CMAPPINGResult> best_result(1000005);
  FILE * fin = fopen(argv[1], "r");
  char cline[MAX_LINE_LENGTH];
  char rname[MAX_LINE_LENGTH];
  //char qname[MAX_LINE_LENGTH];
  // char seq[MAX_LINE_LENGTH];
  //char score[MAX_LINE_LENGTH];
  // char strand, strand2;
  uint32_t start_pos, start_pos2;
  int mismatch, mismatch2;
  cerr << "read ground truth..." << endl;
  uint32_t read = 1;
  uint32_t num_of_paried = 0;
  map<int, int> count_on_mismatch;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    sscanf(cline, "%s %u %u %u %u", rname, &start_pos, &mismatch, &start_pos2,
           &mismatch2);
    best_result[read] = CMAPPINGResult(rname, start_pos, mismatch, start_pos2,
                                       mismatch2);

    if (strcmp(rname, "XXX") != 0 && mismatch != 100) {
      num_of_paried++;
      count_on_mismatch[mismatch + mismatch2]++;
    }
    read++;
  }
  int sum = 0;
  for (map<int, int>::iterator it = count_on_mismatch.begin();
      it != count_on_mismatch.end(); ++it) {
    cout << it->first << " " << it->second << endl;
    sum += it->second;
  }
  cout << "sum = " << sum << endl;
  fclose(fin);

  cerr << "read mapping results..." << endl;
  vector<CMAPPINGResult> res(1000005);
  cerr << argv[3] << endl;
  if (strcmp(argv[2], "bratsw") == 0 || strcmp(argv[2], "-bratsw") == 0) {
    Read_BRATBWResults_pair(argv[3], res);
  } else {
    Read_SAM_Results(argv[3], res, argv[2], true);
  }
  for (size_t i = 0; i < 1000005; ++i) {
    cout << "read: " << i << " (" << res[i].chrom << " " << res[i].start_pos
        << " " << res[i].mismatch << " " << ", " << res[i].start_pos2 << " "
        << res[i].mismatch2 << ")" << endl;
  }

  CompareMappingResultsMR(res, best_result, num_of_paried, count_on_mismatch,
                          sum);

  return 0;
}
