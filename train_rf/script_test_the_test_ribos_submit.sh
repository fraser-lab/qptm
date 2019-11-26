#! /bin/bash

cd /wynton/group/fraser/iyoung/qptm/test_ribo
for tag in `ls`
  do echo "submitting qptxm processing for expt dataset $tag"
  cd $tag
  qsub /wynton/home/fraserlab/iyoung/projects/qptm/script_test_the_test_ribos.sh $tag
  cd /wynton/group/fraser/iyoung/qptm/test_ribo
done


