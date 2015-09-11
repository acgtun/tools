#! /bin/sh
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1:sl230s
#PBS -l mem=35GB
#PBS -m ae
#PBS -q cmb
#PBS -d .

export PATH=$PATH:/home/rcf-40/haifengc/panfs/github/tools/bin

# SRR1532534
./walt -i hg19.dbindex \
       -r ~/panfs/bisulfite_test_data/SRR1532534_1_50000000.fastq \
       -o SRR1532534_single_accu.mr
	   
./walt -i hg19.dbindex \
       -1 ~/panfs/bisulfite_test_data/SRR1532534_1_50000000.fastq \
       -2 ~/panfs/bisulfite_test_data/SRR1532534_2_50000000.fastq \
       -o SRR1532534_paired_accu.mr
	   
accuracy-for-bsmappers ~/panfs/bisulfite_accuracy/SRR1532534/SRR1532534_1_single_benchmark.txt -bsmapper SRR1532534_single_accu.mr
accuracy-for-bsmappers-paired ~/panfs/bisulfite_accuracy/SRR1532534/SRR1532534_pair_benchmark.txt -bsmapper SRR1532534_paired_accu.mr 


# SRR1036985
./walt -i hg19.dbindex \
       -r ~/panfs/bisulfite_test_data/SRR1036985_1_50000000.fastq \
       -o SRR1036985_single_accu.mr
./walt -i hg19.dbindex \
       -1 ~/panfs/bisulfite_test_data/SRR1036985_1_50000000.fastq \
       -2 ~/panfs/bisulfite_test_data/SRR1036985_2_50000000.fastq \
       -o SRR1036985_paired_accu.mr
accuracy-for-bsmappers ~/panfs/bisulfite_accuracy/SRR1036985/SRR1036985_1_single_benchmark.txt -bsmapper SRR1036985_single_accu.mr
accuracy-for-bsmappers-paired ~/panfs/bisulfite_accuracy/SRR1036985/SRR1036985_pair_benchmark.txt -bsmapper SRR1036985_paired_accu.mr

# SRR948855
./walt -i hg19.dbindex \
       -r ~/panfs/bisulfite_test_data/SRR948855_1_50000000.fastq \
       -o SRR948855_single_accu.mr
./walt -i hg19.dbindex \
       -1 ~/panfs/bisulfite_test_data/SRR948855_1_50000000.fastq \
       -2 ~/panfs/bisulfite_test_data/SRR948855_2_50000000.fastq \
       -o SRR948855_paired_accu.mr
	   
accuracy-for-bsmappers ~/panfs/bisulfite_accuracy/SRR948855/SRR948855_1_single_benchmark.txt -bsmapper SRR948855_single_accu.mr
accuracy-for-bsmappers-paired ~/panfs/bisulfite_accuracy/SRR948855/SRR948855_pair_benchmark.txt -bsmapper SRR948855_paired_accu.mr

