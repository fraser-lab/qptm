#! /bin/bash
#$ -j y -cwd -N train_rf -o train.out

time /wynton/home/fraserlab/iyoung/miniconda3/bin/python /wynton/home/fraserlab/iyoung/projects/qptm/train_random_forest_sql.py /wynton/group/fraser/iyoung/qptm/training_combined_results.db

