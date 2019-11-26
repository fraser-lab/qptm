#! /bin/env/bash
#$ -j y -cwd -N prepare -o prep.out

export tag=$1
export iter=$2
export Bmax=$3

cd /wynton/group/fraser/iyoung/qptm/test_ribo/${tag}
for dmin in 2 2.2 2.4 2.6 2.8 3.2 3.8 5
  do for B in `seq 10 10 $Bmax`
    do cd syn_new_dmin${dmin}_B${B}_${iter}
    python3 ~/projects/qptm/prepare_results_for_training_sql.py all_tested_ptms.out synthetic_ptms.out
    cd /wynton/group/fraser/iyoung/qptm/test_ribo/${tag}
  done
done
