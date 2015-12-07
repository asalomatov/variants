#!/bin/bash

cmd=$1
myfile=$2

echo $cmd
echo $myfile
full_script="$cmd $myfile"
echo $full_script

sbatch -J $myfile -N 1 --exclusive -D ./ -o ${myfile}.out -e ${myfile}.err --wrap="$full_script"
