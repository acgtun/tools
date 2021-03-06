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

struct QSUB_OUT {
  QSUB_OUT(const uint32_t& _cputime, const uint32_t& _walttime,
           const double& _memory, const double& _vmem, const string& _file)
      : cputime(_cputime),
        walttime(_walttime),
        memory(_memory),
        vmem(_vmem),
        file(_file) {
  }
  static bool QSUB_OUT_CMP(const QSUB_OUT& a, const QSUB_OUT& b) {
    return a.file < b.file;
  }

  void Output() {
    cout << cputime << "\t" << walttime << "\t" << memory << "\t" << vmem
         << "\t" << file << endl;
  }

  uint32_t cputime;
  uint32_t walttime;
  double memory;
  double vmem;
  string file;
};

int main(int argc, const char *argv[]) {
  string files;
  files = argv[1];
  vector<string> file_names;
  if (isdir(files.c_str())) {
    read_dir(files, file_names);
  } else {
    file_names.push_back(files);
  }

  vector<QSUB_OUT> qsub_out;
  for (uint32_t i = 0; i < file_names.size(); ++i) {
    if (file_names[i].find_first_of('o') == string::npos)
      continue;
    ifstream fin(file_names[i].c_str());
    string line;
    int line_count = 0;
    bool is_pbs_output_file = false;
    while (getline(fin, line)) {
      line_count++;

      if (line_count == 1)
        continue;
      if (line_count == 2) {
        if (line.find("Begin PBS Prologue") != string::npos) {
          is_pbs_output_file = true;
          continue;
        } else {
          break;
        }
      }

      if (line.find("Resources:") == string::npos)
        continue;

      if (is_pbs_output_file == false)
        break;
      uint32_t pos = line.find("cput");
      line = line.substr(pos);
      int hour, min, sec, energy;
      long long mem, vmem;
      sscanf(line.c_str(), "cput=%d:%d:%d,energy_used=%d,mem=%lldkb,vmem=%lld", &hour, &min,
             &sec, &energy, &mem, &vmem);
      uint32_t num_of_seconds = 0;
      num_of_seconds = hour * 3600 + min * 60 + sec;

      ////////////////////////////////////
      uint32_t pwall = line.find("walltime");
      line = line.substr(pwall);
      sscanf(line.c_str(), "walltime=%d:%d:%d", &hour, &min, &sec);
      uint32_t num_of_seconds_wall = 0;
      num_of_seconds_wall = hour * 3600 + min * 60 + sec;
      /////////////////////////////////////////////

      uint32_t p = file_names[i].find_last_of('/');
      if (p != string::npos) {
        file_names[i] = file_names[i].substr(p + 1);
      }
//      qsub_out.push_back(
//          QSUB_OUT((double) num_of_seconds / 3600.00,
//                   (double) num_of_seconds_wall / 3600.00,
//                   (double) mem / (1024.00 * 1024.00),
//                   (double) vmem / (1024.00 * 1024.00), file_names[i]));
      qsub_out.push_back(
          QSUB_OUT(num_of_seconds,
                    num_of_seconds_wall,
                   (double) mem / (1024.00 * 1024.00),
                   (double) vmem / (1024.00 * 1024.00), file_names[i]));
    }
    fin.close();
  }

  sort(qsub_out.begin(), qsub_out.end(), QSUB_OUT::QSUB_OUT_CMP);
  for (uint32_t i = 0; i < qsub_out.size(); ++i) {
    qsub_out[i].Output();
  }

  return 0;
}
