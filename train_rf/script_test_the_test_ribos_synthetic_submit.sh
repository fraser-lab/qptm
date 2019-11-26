#! /bin/bash

cd /wynton/group/fraser/iyoung/qptm/test_ribo
for tag in `ls`
  do cd $tag
  for iter in 2
    do for dmin in 2 2.2 2.4 2.6 2.8 3.2 3.8 5
      do for B in `seq 10 10 90`
        do mkdir syn_new_dmin${dmin}_B${B}_${iter}
        cd syn_new_dmin${dmin}_B${B}_${iter}
        echo "submitting qptxm processing for syn dataset $tag"
        qsub /wynton/home/fraserlab/iyoung/projects/qptm/script_test_the_test_ribos_synthetic.sh $tag $iter $dmin $B
        cd /wynton/group/fraser/iyoung/qptm/test_ribo/${tag}
      done
    done
  done
  cd /wynton/group/fraser/iyoung/qptm/test_ribo
done

