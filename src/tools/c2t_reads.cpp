#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "option.hpp"
#include "smithlab_os.hpp"

using namespace std;

void C2T(const string& org_read, const uint32_t& read_len, string& read) {
  for (uint32_t i = 0; i < read_len; ++i) {
    if ('N' == org_read[i]) {
      read += 'T';
    } else if ('C' == org_read[i]) {
      read += 'T';
    } else {
      read += org_read[i];
    }
  }
}

void A2G(const string& org_read, const uint32_t& read_len, string& read) {
  for (uint32_t i = 0; i < read_len; ++i) {
    if ('N' == org_read[i]) {
      read += 'G';
    } else if ('A' == org_read[i]) {
      read += 'G';
    } else {
      read += org_read[i];
    }
  }
}

void LoadReadsFromFastqFile(const string& reads_file, const bool& AG_WILDCARD) {
  ifstream fin(reads_file.c_str());
  string output_file = reads_file;
  if (AG_WILDCARD) {
    output_file += ".A2G.fq";
  } else {
    output_file += ".C2T.fq";
  }
  ofstream fout(output_file.c_str());

  string line;
  int line_code = 0;
  while (getline(fin, line)) {
    if (line_code == 1) {
      uint32_t read_len = line.size();
      string read;
      if (AG_WILDCARD) {
        A2G(line, read_len, read);
      } else {
        C2T(line, read_len, read);
      }
      fout << read << endl;
    } else {
      fout << line << endl;
    }

    ++line_code;
    if (line_code == 4) {
      line_code = 0;
    }
  }
  fin.close();
  fout.close();
}

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);
  bool AG_WILDCARD = false;
  string reads_file;
  Option::GetOption("-r", reads_file);
  Option::ChkStrExist("-A", AG_WILDCARD);

  LoadReadsFromFastqFile(reads_file, AG_WILDCARD);

  return 0;
}
