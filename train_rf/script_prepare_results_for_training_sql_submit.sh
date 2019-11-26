#! /bin/bash

cd /wynton/group/fraser/iyoung/qptm/test_ribo
for tag in `ls`
  do cd $tag
  echo "prepping results for syn dataset $tag"
  qsub /wynton/home/fraserlab/iyoung/projects/qptm/script_prepare_results_sql.sh $tag 2 90
  cd /wynton/group/fraser/iyoung/qptm/test_ribo
done
