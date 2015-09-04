/*
 * trim-reads.cpp
 *
 *  delete the end of each read, make them shorter
 *
 */

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

#define MAX_LINE_LENGTH 1000

using namespace std;

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);

  string fastq_file;
  string trim_fastq_file;
  uint32_t after_read_len;
  uint32_t num_of_reads;

  Option::GetOption("-f", fastq_file);
  Option::GetOption("-o", trim_fastq_file);
  Option::GetOption("-l", after_read_len, 90);
  Option::GetOption("-n", num_of_reads, 10000000);
  
  FILE * fin = fopen(fastq_file.c_str(), "r");
  ofstream fout(trim_fastq_file.c_str());

  char cline[MAX_LINE_LENGTH];
  string line;
  int line_code = 0;
  uint32_t num = 0;
  while(num < num_of_reads && fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    line = cline;
    switch (line_code) {
      case 0: {
	fout << line << endl;
        break;
      }
      case 1: {
	if(line.size() < after_read_len) {
		std::cout << "read length: " << line.size() << endl;
		exit(EXIT_FAILURE);
	}
        fout << line.substr(0, after_read_len) << endl;
        break;
      }
      case 2: {
	fout << line << endl;
        break;
      }
      case 3: {
	fout << line.substr(0, after_read_len) << endl;
	num++;
        break;
      }
    }
    line_code++; 
    if (line_code == 4) {
      line_code = 0;
    }
  }

  fclose(fin);
  fout.close();

  return 0;
}

