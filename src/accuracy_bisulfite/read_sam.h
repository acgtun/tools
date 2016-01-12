#ifndef TOOLS_READ_SAM_H_
#define TOOLS_READ_SAM_H_

#include <map>
#include <set>
#include <limits>
#include <vector>
#include <cstring>
#include <cstdint>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "option.hpp"
#include "smithlab_os.hpp"

#include <tr1/unordered_map>

using std::tr1::unordered_map;
using namespace std;

const int MAX_LINE_LENGTH = 10000;
const int MAX_INT = std::numeric_limits<int>::max();

struct CMAPPINGResult {
  CMAPPINGResult(string _chrom = "XXX", uint32_t _start_pos = 0,
                 uint32_t _mismatch = 100, uint32_t _start_pos2 = 0,
                 uint32_t _mismatch2 = 100) {
    chrom = _chrom;
    start_pos = _start_pos;
    mismatch = _mismatch;

    start_pos2 = _start_pos2;
    mismatch2 = _mismatch2;
  }

  string chrom;
  uint32_t start_pos;
  uint32_t mismatch;

  uint32_t start_pos2;
  uint32_t mismatch2;
};

struct FLAGLAB {
  bool paired;
  bool paired_mapped;
  bool unmapped;
  bool next_unmapped;
  bool rev;
  bool next_rev;
  bool first;
  bool last;
  bool secondary_align;
};

FLAGLAB GetFLAGLAB(const uint32_t& flag) {
  FLAGLAB flag_lab;
  flag_lab.paired = flag & 0x1;
  flag_lab.paired_mapped = flag & 0x2;
  flag_lab.unmapped = flag & 0x4;
  flag_lab.next_unmapped = flag & 0x8;
  flag_lab.rev = flag & 0x10;
  flag_lab.next_rev = flag & 0x20;
  flag_lab.first = flag & 0x40;
  flag_lab.last = flag & 0x80;
  flag_lab.secondary_align = flag & 0x100;

  return flag_lab;
}

size_t get_mismatch_bismark(const int &mismatch, const string &meth_call_str) {
  /*
   the result of this function might not be accurate, because if a sequencing
   error occurs on a cytosine, then it probably will be reported as a convertion
   */

  int convert_count = 0;
  const char *temp = meth_call_str.substr(5).c_str();
  while (*temp != '\0') {
    if (*temp == 'x' || *temp == 'h' || *temp == 'z')
      ++convert_count;
    ++temp;
  }

  return mismatch - convert_count;
}

void Read_SAM_Results(const char* file_name, vector<CMAPPINGResult>& res,
                      const string& mapper) {
  FILE * fin = fopen(file_name, "r");
  if (!fin) {
    cerr << "cannot open input file " << file_name << endl;
  }

  string qname, rname, cigar, rnext, seq, qual, NM, MD, XM, XR, XG;
  uint32_t flag, mapq, pnext, tlen, pos;
  int mismatch = 100;
  //char strand = '+';

  string qname2, rname2, cigar2, rnext2, seq2, qual2, NM2, MD2, XM2, XR2, XG2;
  uint32_t flag2, mapq2, pnext2, tlen2, pos2 = 0, end_pos2 = 0;
  int mismatch2 = 100;
  //char strand2 = '+';

  char cline[MAX_LINE_LENGTH];
  unsigned int SRRName, SRRName2, line_count = 0;
  while (fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    if (cline[0] == '@')
      continue;
    istringstream iss(cline);
    mismatch = 100;
    mismatch2 = 100;
    pos = 0;
    pos2 = 0;
    if (mapper == "-bsmap" || mapper == "bsmap") {
      iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext
          >> tlen >> seq >> qual >> NM >> MD;
    } else if (mapper == "-bismark" || mapper == "bismark") {
      iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext
          >> tlen >> seq >> qual >> NM >> MD >> XM >> XR >> XG;
    } else if (mapper == "-walt" || mapper == "walt") {
      iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext
          >> tlen >> seq >> qual >> NM;
    } else {
      cout << "the mapper is not correct~!" << endl;
      exit(0);
    }
    sscanf(NM.c_str(), "NM:i:%d", &mismatch);
    if (mapper == "-bismark" || mapper == "bismark") {
      mismatch = get_mismatch_bismark(mismatch, XM);
      if (mismatch2 < 0)
        printf("error~");
    }
    unsigned int read;
    sscanf(qname.c_str(), "SRR%u.%u", &SRRName, &read);
    if (line_count < 1) {
      cout << SRRName << endl;
      cout << cline << endl;
      cout << "mismatch = " << mismatch << endl;
      cout << "pos = " << pos << endl;
    }
    line_count++;
    if (read > 1000000)
    continue;
    FLAGLAB flag_lab = GetFLAGLAB(flag);
    ////////////////////////////////////
    pos2 = 0;
    end_pos2 = 0;
    if (flag_lab.paired_mapped) {
      fgets(cline, MAX_LINE_LENGTH, fin);
      cline[strlen(cline) - 1] = 0;
      if (cline[0] == '@')
        continue;
      istringstream iss(cline);
      //cout << mapper << endl;
      if (mapper == "-bsmap" || mapper == "bsmap") {
        iss >> qname2 >> flag2 >> rname2 >> pos2 >> mapq2 >> cigar2 >> rnext2
            >> pnext2 >> tlen2 >> seq2 >> qual2 >> NM2 >> MD2;
      } else if (mapper == "-bismark" || mapper == "bismark") {
        iss >> qname2 >> flag2 >> rname2 >> pos2 >> mapq2 >> cigar2 >> rnext2
            >> pnext2 >> tlen2 >> seq2 >> qual2 >> NM2 >> MD2 >> XM2 >> XR2
            >> XG2;
      } else if (mapper == "-walt" || mapper == "walt") {
        iss >> qname2 >> flag2 >> rname2 >> pos2 >> mapq2 >> cigar2 >> rnext2
            >> pnext2 >> tlen2 >> seq2 >> qual2 >> NM2;
      } else {
        cout << "the mapper is not correct~!" << endl;
        exit(0);
      }
      sscanf(NM2.c_str(), "NM:i:%d", &mismatch2);
      if (mapper == "-bismark" || mapper == "bismark") {
        mismatch2 = get_mismatch_bismark(mismatch2, XM2);
        if (mismatch2 < 0)
          printf("error~");
      }
      unsigned int read2;
      sscanf(qname2.c_str(), "SRR%u.%u", &SRRName2, &read2);
    }
    res[read] = CMAPPINGResult(rname, pos, mismatch, pos2, mismatch2);
  }

  /* for (int i = 0; i < 1000000; ++i) {
   cout << i << " " << res[i].chrom << " " << res[i].start_pos << " "
   << res[i].mismatch << " " << res[i].start_pos2 << " "
   << res[i].mismatch2 << endl;
   }
   */
  fclose(fin);
}

#endif /* TOOLS_READ_SAM_H_ */
