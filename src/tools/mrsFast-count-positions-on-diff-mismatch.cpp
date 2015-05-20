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
  char chrom[1000];
  unsigned int pos, mismatch, read;
  char strand;
  uint64_t line_count = 0;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    if (line_count % 10000000 == 0) {
      cerr << line_count << endl;
    }
    cline[strlen(cline) - 1] = 0;
    //cout << cline << endl;
    sscanf(cline, "SRR1171540.%u %s %u %c %u", &read, chrom, &pos,
           &strand, &mismatch);
    //cout << read << " " << mismatch << endl;
    count[read][mismatch]++;
    line_count++;
  }
  fclose(fin);

  string file = argv[1];
  file += ".count.txt";
  FILE * fout = fopen(file.c_str(), "w");
  for (uint32_t i = 0; i <= 1000000; ++i) {
    fprintf(fout, "%u:", i);
    for (uint32_t j = 0; j <= 6; ++j) {
      fprintf(fout, " %u", count[i][j]);
    }
    fprintf(fout, "\n");
  }

  fclose(fout);
}
