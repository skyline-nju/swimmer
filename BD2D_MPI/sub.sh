#!/bin/bash

L=1500
phi=0.04
n_step=5000000
h=0.005
k=5

for alpha in 0.005 0.006 0.007 0.008 0.009 0.01 0.013 0.016 0.02 0.025 0.03
do
    f=${alpha}_${phi}_${k}
    sub -N $f ./main.out -L ${L} --phi ${phi} -h ${h} --alpha ${alpha} -k ${k} -o data --profile --log_dt 50000 --profile_dt 1000 --snap_dt 10000 --n_step ${n_step}
    sleep 1
done
