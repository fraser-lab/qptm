#! /bin/bash

cd /wynton/group/fraser/iyoung/qptm/test_ribo
for tag in `ls`
  do cd $tag
  echo "prepping results for syn dataset $tag"
  qsub /wynton/home/fraserlab/iyoung/projects/qptm/script_prepare_results_sql.sh $tag 3 10
  cd /wynton/group/fraser/iyoung/qptm/test_ribo
done

cd /wynton/group/fraser/iyoung/qptm/test_ribo
for tag in `ls`
  do cd $tag
  echo "prepping results for expt dataset $tag"
  python3 ~/projects/qptm/prepare_results_for_training_sql.py all_tested_ptms.out modeled_ptms.out
  cd /wynton/group/fraser/iyoung/qptm/test_ribo
done

cd /wynton/group/fraser/iyoung/qptm/final_ribo
for iter in 2
  do for dmin in 2 2.2 2.4 2.6 2.8 3.2 3.8 5
    do for B in 10
      do cd syn_new_dmin${dmin}_B${B}_${iter}
      echo "prepping results for syn dataset for Stojkovic ribo"
      python3 ~/projects/qptm/prepare_results_for_training_sql.py all_tested_ptms.out synthetic_ptms.out
      cd /wynton/group/fraser/iyoung/qptm/final_ribo
    done
  done
done

cd /wynton/group/fraser/iyoung/qptm/final_ribo
echo "prepping results for expt Stojkovic ribo"
python3 ~/projects/qptm/prepare_results_for_training_sql.py all_tested_ptms.out modeled_ptms.out

