#! /bin/bash

python3 ~/projects/qptm/combine_results_for_training_sql.py `find /wynton/group/fraser/iyoung/qptm/test_ribo/*/syn*_2/ -name "results.db"`
mv combined_results.db training_combined_results.db

