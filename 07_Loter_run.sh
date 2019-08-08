#!/bin/bash

for j in 10 20 30 40 50 75 100
do
echo $j
~/software/Local_Ancestry/Loter.sh \
--ref ~/data/sim_${j}/intermediate/sim_${j}.ref.bcf.gz \
--query ~/data/sim_${j}/intermediate/sim_${j}.query.vcf \
--pop ~/data/sim_${j}/intermediate/sim_${j}.ref.map \
--out ~/data/sim_${j}/loter_benchmarking_sim_${j}
done
