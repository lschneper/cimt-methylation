#!/bin/bash
#cd /drives/drive1/ff_chaya/Lisa/ewas_results/y22/y22sensitivity
#for i in {2..7}
#do

Rscript /drives/drive1/ff_chaya/Lisa/April_2026/bin/Analysis/sensitivity/y22_epic_sex_ewas_models_042026.R 1 100000
Rscript /drives/drive1/ff_chaya/Lisa/April_2026/bin/Analysis/sensitivity/y22_epic_sex_ewas_models_042026.R 100001 200000
Rscript /drives/drive1/ff_chaya/Lisa/April_2026/bin/Analysis/sensitivity/y22_epic_sex_ewas_models_042026.R 200001 300000
Rscript /drives/drive1/ff_chaya/Lisa/April_2026/bin/Analysis/sensitivity/y22_epic_sex_ewas_models_042026.R 300001 400000
Rscript /drives/drive1/ff_chaya/Lisa/April_2026/bin/Analysis/sensitivity/y22_epic_sex_ewas_models_042026.R 400001 500000
Rscript /drives/drive1/ff_chaya/Lisa/April_2026/bin/Analysis/sensitivity/y22_epic_sex_ewas_models_042026.R 500001 600000
Rscript /drives/drive1/ff_chaya/Lisa/April_2026/bin/Analysis/sensitivity/y22_epic_sex_ewas_models_042026.R 600001 700000
Rscript /drives/drive1/ff_chaya/Lisa/April_2026/bin/Analysis/sensitivity/y22_epic_sex_ewas_models_042026.R 700001 800000
Rscript /drives/drive1/ff_chaya/Lisa/April_2026/bin/Analysis/sensitivity/y22_epic_sex_ewas_models_042026.R 800001 900000
#done

#cat y22_epic_start_cpg_100000_model${i}.ewasoutput.txt > y22_epic_model${i}.ewasoutput.txt
#awk '{if(NR>1){print $0}}' y22_epic_start_cpg_200000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
#awk '{if(NR>1){print $0}}' y22_epic_start_cpg_300000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
#awk '{if(NR>1){print $0}}' y22_epic_start_cpg_400000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
#awk '{if(NR>1){print $0}}' y22_epic_start_cpg_500000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
#awk '{if(NR>1){print $0}}' y22_epic_start_cpg_600000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
#awk '{if(NR>1){print $0}}' y22_epic_start_cpg_700000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
#awk '{if(NR>1){print $0}}' y22_epic_start_cpg_800000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
#awk '{if(NR>1){print $0}}' y22_epic_start_cpg_900000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
#done
