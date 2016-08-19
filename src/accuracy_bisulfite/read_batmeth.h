#include "read_sam.h"

void Read_BatMethResults(const char* file_name, vector<CMAPPINGResult>& res) {
  FILE * fin = fopen(file_name, "r");
  if (!fin) {
    cerr << "cannot open input file " << file_name << endl;
  }

  string qname, rname, strand;
  uint32_t pos, read;
  int mismatch = 100;

  char cline[MAX_LINE_LENGTH];
  unsigned int SRRName, SRRName2, line_count = 0, ore = 0;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    if (cline[0] == '@' && strlen(cline) == 1)
      continue;
    if (cline[0] == '@') {
      istringstream iss(cline);
      cout << cline << endl;
      sscanf(cline, "@SRR%u.%u", &SRRName, &read);
      cout << read << endl;
      if (read > 1000000)
        break;
      map<int, vector<CMAPPINGResult> > results;
      while (fgets(cline, MAX_LINE_LENGTH, fin)) {
        cline[strlen(cline) - 1] = 0;
        if (cline[0] == '@') {
          cout << "---------------------------" << endl;
          for (map<int, vector<CMAPPINGResult> >::iterator it = results.begin();
                 it != results.end(); ++it) {
            cout << it->first << " " << it->second.size() << endl;
            for (size_t i = 0;i < it->second.size();++i) {
              cout << it->second[i].chrom << " " << it->second[i].start_pos << " " << it->second[i].mismatch << endl;
            }
           }
          for (map<int, vector<CMAPPINGResult> >::iterator it = results.begin();
              it != results.end(); ++it) {
            if (it->second.size() > 1)
              break;
            if (it->second.size() == 1) {
              res[read] = it->second[0];
              break;
            }
          }
          break;
        }

        //////
        istringstream iss(cline);
        iss >> ore >> rname >> strand >> pos >> mismatch;
        results[mismatch].push_back(
            CMAPPINGResult(rname, pos + 1, mismatch, 0, 100));
      }
    }
  }

  fclose(fin);
}
