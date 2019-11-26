#! /bin/bash
#$ -j n -cwd -N syn -o log.out -e log.err
. /wynton/home/fraserlab/iyoung/apps/phenix/phenix-1.17.1-3660/build/setpaths.sh

export iter=$1
export dmin=$2
export B=$3

export model="4ybb-pdb-bundle2_ali-coot-40_PSU_as_U.pdb"
export map="50S_Ecoli_big_2p15A-noBfau.mrc"
echo "testing dmin=${dmin}, B=${B}, iter ${iter}..."
cp ../${model} .
cp ../${map} .
echo "executing libtbx.python /wynton/home/fraserlab/iyoung/projects/qptm/qptxm.py map_file=${map} model_file=${model} d_min=${dmin} set_b_factor=${B} plot=False synthetic_data=True"
libtbx.python /wynton/home/fraserlab/iyoung/projects/qptm/qptxm.py map_file=${map} model_file=${model} d_min=${dmin} set_b_factor=${B} plot=False synthetic_data=True
rm /wynton/group/fraser/iyoung/qptm/final_ribo/syn_new_dmin${dmin}_B${B}_${iter}/$model 
rm /wynton/group/fraser/iyoung/qptm/final_ribo/syn_new_dmin${dmin}_B${B}_${iter}/$map
rm *.ccp4
touch "done"

