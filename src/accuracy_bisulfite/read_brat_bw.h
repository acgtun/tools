#include "read_sam.h"

void Read_BRATBWResults(const char* file_name, vector<CMAPPINGResult>& res) {
  FILE * fin = fopen(file_name, "r");
  if (!fin) {
    cerr << "cannot open input file " << file_name << endl;
  }

  string qname, rname, qstr, strand;
  uint32_t pos, read, pos2;
  int mismatch = 100;

  char cline[MAX_LINE_LENGTH];
  unsigned int SRRName, SRRName2, line_count = 0, ore = 0;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    cout << cline << endl;
    istringstream iss(cline);
    iss >> read >> qstr >> rname >> strand >> pos >> mismatch >> pos2;
    read = read + 1;
    cout << "haha: " << read << " " << strand << " " << pos << " " << mismatch
        << endl;
    if(read > 1000000) continue;
    res[read] = CMAPPINGResult(rname, pos + 1, mismatch, pos2 + 1, 100);
  }

  fclose(fin);
}

void Read_BRATBWResults_pair(const char* file_name, vector<CMAPPINGResult>& res) {
  FILE * fin = fopen(file_name, "r");
  if (!fin) {
    cerr << "cannot open input file " << file_name << endl;
  }

  string qname, rname, qstr, qstr2, strand;
  uint32_t pos, read, pos2, pos_f, pos_f2, mismatch1, mismatch2;
  int mismatch = 100;

  char cline[MAX_LINE_LENGTH];
  unsigned int SRRName, SRRName2, line_count = 0, ore = 0;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    cout << cline << endl;
    istringstream iss(cline);
    iss >> read >> qstr  >> qstr2 >> rname >> strand >> pos_f >> pos_f2 >> mismatch1 >> mismatch2 >> pos >> pos2;
    read = read + 1;
    cout << "haha: " << read << " " << strand << " " << pos << " " << mismatch1
        << endl;
    if(read > 1000000) continue;
    res[read] = CMAPPINGResult(rname, pos_f + 1, mismatch1, pos + 1, mismatch2);
  }

  fclose(fin);
}
