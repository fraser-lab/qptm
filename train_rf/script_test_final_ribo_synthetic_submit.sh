#! /bin/bash

cd /wynton/group/fraser/iyoung/qptm/final_ribo
export iter=2
for dmin in 2 2.2 2.4 2.6 2.8 3.2 3.8 5
  do for B in `seq 10 10 90`
    do mkdir syn_new_dmin${dmin}_B${B}_${iter}
    cd syn_new_dmin${dmin}_B${B}_${iter}
    echo "submitting qptxm processing for synthetic Stojkovic ribosome"
    qsub /wynton/home/fraserlab/iyoung/projects/qptm/script_test_final_ribo_synthetic.sh $iter $dmin $B
    cd /wynton/group/fraser/iyoung/qptm/final_ribo
  done
done

