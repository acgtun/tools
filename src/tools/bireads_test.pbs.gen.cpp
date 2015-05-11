#include <string>
#include <fstream>
#include <iostream>

using namespace std;

string folder[] = { "SRX472622", "SRX666252", "SRP041996" };
int read_files_size[3] = {3, 6, 3};
int maxL[3] = {29, 29, 36};
int N[] = {1000000, 5000000};
string reads_files[3][6] = {
    {"SRR1171540", "SRR1171541", "SRR1171542"},
    {"SRR1532534", "SRR1532535", "SRR1532536", "SRR1532537", "SRR1532538", "SRR1532539"},
    {"SRR1287225", "SRR1287226","SRR1287227"}
};

int main(){
  for(int i = 0;i < 3;++i) {
    for(int j = 0;j < read_files_size[i];++j){
      for(int L = 25;L <= maxL[i];++L){
        for(int n = 0;n < 2;++n){
          for(int p = 1;p <= 3;++p){
            char file[100];
            sprintf(file, "%s_%d_%d_%d_bsmapper.pbs", reads_files[i][j].c_str(), p, L, N[n]);
            ofstream fout(file);
            fout << "#! /bin/sh" << endl;
            fout << "#PBS -l walltime=500:00:00" << endl;
            fout << "#PBS -l nodes=1:ppn=1:sl230s" << endl;
            fout << "#PBS -l mem=25GB" << endl;
            fout << "#PBS -q cmb" << endl;
            fout << "#PBS -d ." << endl;
            fout << endl;

            fout << "./bsmapper \\" << endl;
            fout << "      -i ~/panfs/github/smithlab_cpp/hg19_index/hg19.dbindex \\" << endl;

            if(p == 3){
              fout << "      -1 ~/panfs/bireads/" << folder[i] << "/" << reads_files[i][j] << "_1.fastq \\" << endl;
              fout << "      -2 ~/panfs/bireads/" << folder[i] << "/" << reads_files[i][j] << "_2.fastq \\" << endl;
            } else {
              fout << "      -r ~/panfs/bireads/" << folder[i] << "/" << reads_files[i][j] << "_" << p << ".fastq \\" << endl;
            }

            fout << "      -o " << reads_files[i][j] << "_" << p << "_L" << L << "_" << N[n] << ".out \\" << endl;
            fout << "      -l " << L  << " \\"<< endl;
            fout << "      -N " << N[n] << endl;

            fout << endl;

            fout.close();
          }
        }
      }
    }
  }
}

