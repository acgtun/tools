#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "smithlab_os.hpp"

using namespace std;

void single_end_run(const string& fq_file) {
  char file[100];
  sprintf(file, "%s_single.pbs", fq_file.c_str());
  ofstream fout(file);
  fout << "#! /bin/sh" << endl;
  fout << "#PBS -l walltime=200:00:00" << endl;
  fout << "#PBS -l nodes=1:ppn=1:sl230s" << endl;
  fout << "#PBS -l mem=18GB" << endl;
  fout << "#PBS -q cmb" << endl;
  fout << "#PBS -d ." << endl;
  fout << endl;

  fout << "/home/rcf-40/haifengc/panfs/bin/walt \\" << endl;
  fout << "      -i /home/rcf-40/haifengc/panfs/hg19_walt_index_random/hg19.dbindex \\"
       << endl;
  fout << "      -r " << fq_file << " \\" << endl;
  fout << "      -o " << fq_file << "_walt_single_out.mr" << endl;
  fout << endl;
  fout.close();
}

void paired_end_run(const string& fq_file1, const string& fq_file2) {
  char file[100];
  sprintf(file, "%s_pair.pbs", fq_file1.c_str());
  ofstream fout(file);
  fout << "#! /bin/sh" << endl;
  fout << "#PBS -l walltime=200:00:00" << endl;
  fout << "#PBS -l nodes=1:ppn=1:sl230s" << endl;
  fout << "#PBS -l mem=28GB" << endl;
  fout << "#PBS -q cmb" << endl;
  fout << "#PBS -d ." << endl;
  fout << endl;

  fout << "/home/rcf-40/haifengc/panfs/bin/walt \\" << endl;
  fout << "      -i /home/rcf-40/haifengc/panfs/hg19_walt_index_random/hg19.dbindex \\"
       << endl;
  fout << "      -1 " << fq_file1 << " \\" << endl;
  fout << "      -2 " << fq_file2 << " \\" << endl;
  fout << "      -o " << fq_file1 << "_walt_pair_out.mr" << endl;
  fout << endl;
  fout.close();
}

bool check_paired_end(const string& fq_file1, const string& fq_file2) {
  size_t p1 = fq_file1.find("_1.f");
  size_t p2 = fq_file2.find("_2.f");
  if (p1 != string::npos && p2 != string::npos
      && fq_file1.substr(0, p1) == fq_file2.substr(0, p2))
    return true;

  return false;
}

int main(int argc, const char *argv[]) {
  string files;
  files = argv[1];
  vector<string> file_names;
  if (isdir(files.c_str())) {
    read_dir(files, file_names);
  } else {
    file_names.push_back(files);
  }

  sort(file_names.begin(), file_names.end());

  for (uint32_t i = 0; i < file_names.size(); ++i) {
    cout << file_names[i] << endl;
    if (i + 1 >= file_names.size()
        && (is_valid_filename(file_names[i], "fastq")
            || is_valid_filename(file_names[i], "fq"))) {
      single_end_run(file_names[i]);
    } else if (i + 1 < file_names.size()) {
      if (check_paired_end(file_names[i], file_names[i + 1])) {
        paired_end_run(file_names[i], file_names[i + 1]);
        i++;
      } else if (is_valid_filename(file_names[i], "fastq")
          || is_valid_filename(file_names[i], "fq")) {
        single_end_run(file_names[i]);
      }
    }
  }

  return 0;
}
