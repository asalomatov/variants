#!/bin/bash


child_id=$1 # trio002.p1_641945
cfg_yml=$2 # /mnt/xfs1/home/asalomatov/projects/variants/variants/cfg_spark.yml
model=$3 # model_final_0516/GBM_lvl9_stdFalse_cut0.5_splt0.99_5000__0516_2000_1_0.075.pkl
prob_cutoff=$4
echo $@

echo "running /mnt/xfs1/home/asalomatov/projects/variants/variants/work/extractFeatures.py $child_id $cfg_yml"
/mnt/xfs1/home/asalomatov/miniconda2/bin/python \
    /mnt/xfs1/home/asalomatov/projects/variants/variants/work/extractFeatures.py $child_id $cfg_yml

echo "running /mnt/xfs1/home/asalomatov/projects/variants/variants/test.py $child_id $cfg_yml $model 0 0 $prob_cutoff"
/mnt/xfs1/home/asalomatov/miniconda2/bin/python \
    /mnt/xfs1/home/asalomatov/projects/variants/variants/test.py \
    $child_id $cfg_yml $model 0 0 $prob_cutoff
