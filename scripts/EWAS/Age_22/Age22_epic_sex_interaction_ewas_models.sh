#!/bin/bash
cd <path to Age22_epic_sex_interaction_ewas_model.R>


Rscript Age22_epic_sex_interaction_ewas_model.R 1 100000
Rscript Age22_epic_sex_interaction_ewas_model.R 100001 200000
Rscript Age22_epic_sex_interaction_ewas_model.R 200001 300000
Rscript Age22_epic_sex_interaction_ewas_model.R 300001 400000
Rscript Age22_epic_sex_interaction_ewas_model.R 400001 500000
Rscript Age22_epic_sex_interaction_ewas_model.R 500001 600000
Rscript Age22_epic_sex_interaction_ewas_model.R 600001 700000
Rscript Age22_epic_sex_interaction_ewas_model.R 700001 800000
Rscript Age22_epic_sex_interaction_ewas_model.R 800001 900000

i=9
cat y22_epic_start_cpg_100000_model${i}.ewasoutput.txt > y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_200000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_300000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_400000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_500000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_600000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_700000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_800000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_900000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
