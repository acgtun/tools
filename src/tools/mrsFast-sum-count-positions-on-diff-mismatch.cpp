#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

const int MAX_LINE_LENGTH = 1000;
int main(int argc, const char **argv) {
  FILE * fin = fopen(argv[1], "r");
  vector<vector<uint32_t> > count(1000005, vector<uint32_t>(7, 0));
  char cline[MAX_LINE_LENGTH];
  unsigned int read, m0, m1, m2, m3, m4, m5, m6;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    sscanf(cline, "%u: %u %u %u %u %u %u %u", &read, &m0, &m1, &m2, &m3, &m4,
           &m5, &m6);
    count[read][0] = m0;
    count[read][1] = m1;
    count[read][2] = m2;
    count[read][3] = m3;
    count[read][4] = m4;
    count[read][5] = m5;
    count[read][6] = m6;
  }
  fclose(fin);
  /////////////////////////////////////////////
  fin = fopen(argv[2], "r");
  vector<vector<uint32_t> > count2(1000005, vector<uint32_t>(7, 0));
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    sscanf(cline, "%u: %u %u %u %u %u %u %u", &read, &m0, &m1, &m2, &m3, &m4,
           &m5, &m6);
    count2[read][0] = m0;
    count2[read][1] = m1;
    count2[read][2] = m2;
    count2[read][3] = m3;
    count2[read][4] = m4;
    count2[read][5] = m5;
    count2[read][6] = m6;
  }
  fclose(fin);

  string file = argv[1];
  file += ".sum.txt";
  FILE * fout = fopen(file.c_str(), "w");
  for (uint32_t i = 1; i <= 1000000; ++i) {
    fprintf(fout, "SRR1171540.%u", i);
    for (uint32_t j = 0; j <= 6; ++j) {
      fprintf(fout, " %u", count[i][j] + count2[i][j]);
    }
    fprintf(fout, "\n");
  }

  fclose(fout);
}
