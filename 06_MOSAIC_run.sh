#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00

for pop in sim_10 sim_20 sim_30 sim_40 sim_50 sim_75 sim_100
do

if [ ${pop} == "EVEN" ] || [ ${pop} == "BRYC" ]
then
	nanc=3
else
	nanc=2
fi

if [ ${pop} == "EVEN" ] || [ ${pop} == "BRYC" ] || [ ${pop} == "PUR" ] || [ ${pop} == "MXL" ]
then
	ngen=11
else
	ngen=7
fi
mkdir -p ~/software/Local_Ancestry/MOSAIC_benchmarking/${pop}

~/software/Local_Ancestry/06_MOSAIC.sh \
--ancestries ${nanc} \
--query ~/data/${pop}/intermediate/${pop}.query.vcf \
--ref ~/data/${pop}/intermediate/${pop}.ref.bcf.gz \
--geneticmap ~/data/ancestry_validation/cohort1/thinned/ASW_thinning_thinned_genetic_map_intersection.txt \
--pop ~/data/${pop}/intermediate/${pop}.ref.map \
--generations ${ngen} \
--outdir ~/software/Local_Ancestry/MOSAIC_benchmarking/${pop} \
--prefix ${pop}
done
#/usr/bin/time  ~/enet_scripts/prune_testing/pop_train_combos_unpruned.R ${chr} ${pop}
