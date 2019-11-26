#! /bin/bash

export db="results.db"
export rf="/wynton/group/fraser/iyoung/qptm/trained_rf.pkl"

cd /wynton/group/fraser/iyoung/qptm/test_ribo
for tag in `ls`
  do cd $tag
  qsub /wynton/home/fraserlab/iyoung/projects/qptm/script_plot_results_sql.sh $tag
  cd /wynton/group/fraser/iyoung/qptm/test_ribo
done

cd /wynton/group/fraser/iyoung/qptm/final_ribo/
export title="\"Stojkovic dataset\n2.2 Å resolution\""
echo "plotting results for expt Stojkovic ribo" > plot.out
/wynton/home/fraserlab/iyoung/miniconda3/bin/python /wynton/home/fraserlab/iyoung/projects/qptm/test_random_forest_sql.py $rf $db $title >> plot.out

cd /wynton/group/fraser/iyoung/qptm/final_ribo/
qsub /wynton/home/fraserlab/iyoung/projects/qptm/script_plot_final_ribo_results_sql.sh

