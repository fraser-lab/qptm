#! /bin/bash

find /wynton/group/fraser/iyoung/qptm/test_ribo/*/syn*_3/ -name "results.db" > test_synthetic
find /wynton/group/fraser/iyoung/qptm/test_ribo/ -name "results.db" | grep -v syn > test_experimental
find /wynton/group/fraser/iyoung/qptm/final_ribo/ -name "results.db" > test_final
python3 ~/projects/qptm/combine_results_for_training_sql.py $(cat test_synthetic) $(cat test_experimental) $(cat test_final)
mv combined_results.db testing_combined_results.db

