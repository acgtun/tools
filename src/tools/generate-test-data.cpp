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

int main(int argc, const char *argv[]) {
  InitProgram(argc, argv);
  string database, output_file;
  int query_length, num_of_query;
  bool fix_length, gap_exist;
  Option::GetOption("-d", database);
  Option::GetOption("-l", query_length, 100);

  /* set all queries to have the same length*/
  Option::ChkStrExist("-f", fix_length);

 /* will random generate gaps in the queries*/ 
  Option::ChkStrExist("-g", gap_exist);

  Option::GetOption("-n", num_of_query, 10000);
  Option::GetOption("-o", output_file);

  ifstream fin(database.c_str());

  vector<string> chrom_or_protein_names;
  vector<string> chrom_or_protein_seqs;
  read_fasta_file(database.c_str(), chrom_or_protein_names,
                  chrom_or_protein_seqs);

  int num_of_chrom_or_protein = chrom_or_protein_names.size();
  srand(time(NULL));
  ofstream fout(output_file.c_str());
  for (int i = 0; i < num_of_query; ++i) {
    uint32_t id, start_pos, len;
    while (1) {
      id = rand() % num_of_chrom_or_protein;
      start_pos = rand() % 1000;
      if (fix_length)
        len = query_length;
      else
        len = rand() % 1000 + 6;
      if (start_pos + len < chrom_or_protein_seqs[id].size())
        break;
    }

    fout << ">test" << i + 1 << " " << chrom_or_protein_names[id] << " "
        << start_pos << " " << len << endl;
    for (uint32_t k = start_pos, l = 0; l < len; k++, l++) {
      if (l % 80 == 0 && l != 0) {
        fout << endl;
      }
      int r = rand() % 250;
      if (r == 0 && gap_exist)
        continue;
      fout << chrom_or_protein_seqs[id][k];

    }
    fout << endl;
  }
  return 0;
}

