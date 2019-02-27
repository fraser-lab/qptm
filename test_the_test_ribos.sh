#! /bin/bash

. ~/source_xfel_dev.sh

cd /Users/iris/projects/qptm/test_ribo
for tag in `ls`
  do echo "testing ribosome model ${tag}..."
  cd $tag
  export model=`ls | grep cif$ | tail -1`
  export map=`ls | grep map$ | tail -1`
  export command="libtbx.python /Users/iris/projects/qptm/qptm.py map_file=${map} model_file=${model}"
  echo "executing $command"
  $command
  echo "press return to open the notes for this test"
  read response
  vim test.out
  #echo "add any notes, terminating in $"
  #read -d "$" notes
  #echo "\n" >> test.out
  #echo $notes >> test.out
  #echo
  cd /Users/iris/projects/qptm/test_ribo
done


