#! /bin/sh
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1:sl230s
#PBS -l mem=10GB
#PBS -m ae
#PBS -q cmb
#PBS -d .

export PATH=$PATH:/home/rcf-40/haifengc/panfs/github/smithlabcode/methpipe/bin


mkdir mr
mv *.mr mr

# merge reads mapping results
cat ./mr/*.mr > merge.mr

# sort mapping results
LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o merge.mr.sorted_start merge.mr

# remove duplicates
duplicate-remover -S merge.mr_dremove_stat.txt -o merge.mr.dremove merge.mr.sorted_start

# bsrate
bsrate -c ~/panfs/hg19 -o merge.bsrate merge.mr.dremove

# methcounts
methcounts -c ~/panfs/hg19 -o merge.mr.meth merge.mr.dremove

#levels
levels -o merge.levels merge.mr.meth
