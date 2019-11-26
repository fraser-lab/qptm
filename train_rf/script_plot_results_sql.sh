#! /bin/env/bash
#$ -j y -cwd -N plot -o plot_log.out

export tag=$1
export iter=3
export B=10
export db="results.db"
export rf="/wynton/group/fraser/iyoung/qptm/trained_rf.pkl"

cd /wynton/group/fraser/iyoung/qptm/test_ribo/${tag}
for dmin in 2 2.2 2.4 2.6 2.8 3.2 3.8 5
  do cd syn_new_dmin${dmin}_B${B}_${iter}
  export title="\"synthetic dataset ${tag}, dmin=${dmin} Å, B=${B}\""
  echo "plotting results for syn dataset $tag" > plot.out
  /wynton/home/fraserlab/iyoung/miniconda3/bin/python /wynton/home/fraserlab/iyoung/projects/qptm/test_random_forest_sql.py $rf $db $title >> plot.out
  cd /wynton/group/fraser/iyoung/qptm/test_ribo/${tag}
done

cd /wynton/group/fraser/iyoung/qptm/test_ribo/${tag}
export dmin=$(cat dmin)
export title="\"experimental dataset ${tag}, dmin=${dmin} Å\""
echo "plotting results for expt dataset $tag" > plot.out
/wynton/home/fraserlab/iyoung/miniconda3/bin/python /wynton/home/fraserlab/iyoung/projects/qptm/test_random_forest_sql.py $rf $db $title >> plot.out

