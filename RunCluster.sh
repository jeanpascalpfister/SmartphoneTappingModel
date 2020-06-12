#!/bin/bash
module load matlab

nsim=(11 11 1 1 12 12) # number of simulations per Simulation serie
for serie in {1..6}
do
ns=${nsim[$serie-1]}
list=($(seq 1 1 $ns))
for s in ${list[@]}
do
for p in {1..84}
do
        echo "Queuing up serie=${serie}, sim = ${s}, subject=${p}"
        # Submit the job
bsub -W "24:00" -R "rusage[mem=1024]" -N -o log2/logSER${serie}_SUB${p}_SIM${s}.txt -J "SER=${serie}, SUB=${p}, SIM=${s}" matlab -singleCompThread -nosplash -nodesktop -r "ClusterFitSubject(${p},${serie},${s},'DataSet.mat','fit2/SER${serie}_SUB${p}_SIM${s}.mat'),exit"
done
done
done
wait
echo "Done"
