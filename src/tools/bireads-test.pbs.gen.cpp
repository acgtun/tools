#include <string>
#include <fstream>
#include <iostream>

using namespace std;

string folder[] = { "SRX472622", "SRX666252", "SRP041996" };
int read_files_size[3] = { 3, 6, 3 };
string suffix[] = { "mr", "sam" };
string reads_files[3][6] = {
    { "SRR1171540", "SRR1171541", "SRR1171542" },
    { "SRR1532534", "SRR1532535", "SRR1532536", "SRR1532537", "SRR1532538", "SRR1532539" },
    { "SRR1287225", "SRR1287226", "SRR1287227" }
};

void bismark() {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < read_files_size[i]; ++j) {
      char file[100];
      sprintf(file, "%s_1_bismark.pbs", reads_files[i][j].c_str());
      ofstream fout(file);
      fout << "#! /bin/sh" << endl;
      fout << "#PBS -l walltime=500:00:00" << endl;
      fout << "#PBS -l nodes=1:ppn=1:sl230s" << endl;
      fout << "#PBS -l mem=15GB" << endl;
      fout << "#PBS -q cmb" << endl;
      fout << "#PBS -d ." << endl;
      fout << endl;

      fout
          << "./bismark -N 1 -L 32 --path_to_bowtie ~/panfs/bisulfite_mapping_program/bismark_v0.13.1/bowtie2-2.2.4/ --bowtie2 ~/panfs/bsm_test/bismark/hg19 ~/panfs/bireads/"
          << folder[i] << "/" << reads_files[i][j] << "_1.fastq" << endl;
      fout.close();
    }
  }
}

void bsmap() {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < read_files_size[i]; ++j) {
      char file[100];
      sprintf(file, "%s_1_bsmap.pbs", reads_files[i][j].c_str());
      ofstream fout(file);
      fout << "#! /bin/sh" << endl;
      fout << "#PBS -l walltime=500:00:00" << endl;
      fout << "#PBS -l nodes=1:ppn=1:sl230s" << endl;
      fout << "#PBS -l mem=10GB" << endl;
      fout << "#PBS -q cmb" << endl;
      fout << "#PBS -d ." << endl;
      fout << endl;

      fout << "./bsmap \\" << endl;
      fout << "   -a ~/panfs/bireads/" << folder[i] << "/" << reads_files[i][j]
           << "_1.fastq \\" << endl;
      fout << "   -d ~/panfs/hg19_all.fa \\" << endl;
      fout << "   -v 6 -p 1 -r 0 \\" << endl;
      fout << "   -o bsmap_" << reads_files[i][j] << ".sam" << endl;
      fout.close();
    }
  }
}

int main() {
  //bismark();
  //bsmap();
  //return 0;
  //single
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < read_files_size[i]; ++j) {
      for (int s = 0; s < 2; ++s) {
        char file[100];
        sprintf(file, "%s_single_walt_%s.pbs", reads_files[i][j].c_str(), suffix[s].c_str());
        ofstream fout(file);
        fout << "#! /bin/sh" << endl;
        fout << "#PBS -l walltime=500:00:00" << endl;
        fout << "#PBS -l nodes=1:ppn=1:sl230s" << endl;
        fout << "#PBS -l mem=18GB" << endl;
        fout << "#PBS -q cmb" << endl;
        fout << "#PBS -d ." << endl;
        fout << endl;

        fout << "./walt \\" << endl;
        fout << "      -i hg19.dbindex \\" << endl;
        fout << "      -r ~/panfs/bireads/" << folder[i] << "/" << reads_files[i][j] << "_1.fastq \\" << endl;
        fout << "      -o " << reads_files[i][j] << "_single_1." << suffix[s] << endl;
        fout << endl;

        fout.close();
      }
    }
  }

  // pair
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < read_files_size[i]; ++j) {
      for (int s = 0; s < 2; ++s) {
        char file[100];
        sprintf(file, "%s_pair_walt_%s.pbs", reads_files[i][j].c_str(), suffix[s].c_str());
        ofstream fout(file);
        fout << "#! /bin/sh" << endl;
        fout << "#PBS -l walltime=500:00:00" << endl;
        fout << "#PBS -l nodes=1:ppn=1:sl230s" << endl;
        fout << "#PBS -l mem=28GB" << endl;
        fout << "#PBS -q cmb" << endl;
        fout << "#PBS -d ." << endl;
        fout << endl;

        fout << "./walt \\" << endl;
        fout << "      -i hg19.dbindex \\" << endl;
        fout << "      -1 ~/panfs/bireads/" << folder[i] << "/" << reads_files[i][j] << "_1.fastq \\" << endl;
        fout << "      -2 ~/panfs/bireads/" << folder[i] << "/" << reads_files[i][j] << "_2.fastq \\" << endl;
        fout << "      -o " << reads_files[i][j] << "_pair." << suffix[s] << endl;
        fout << endl;

        fout.close();
      }
    }
  }
  return 0;
}
