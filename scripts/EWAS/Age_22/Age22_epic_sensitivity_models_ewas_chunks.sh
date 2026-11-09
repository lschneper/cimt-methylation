#!/bin/bash
cd <path to Rscript Age22_epic_sensitivity_models_ewas_chunks.R>
for i in {1..7} {10..13};
do
	Rscript Age22_epic_sensitivity_models_ewas_chunks.R 1 100000 ${i}
	Rscript Age22_epic_sensitivity_models_ewas_chunks.R 100001 200000 ${i}
	Rscript Age22_epic_sensitivity_models_ewas_chunks.R 200001 300000 ${i}
	Rscript Age22_epic_sensitivity_models_ewas_chunks.R 300001 400000 ${i}
	Rscript Age22_epic_sensitivity_models_ewas_chunks.R 400001 500000 ${i}
	Rscript Age22_epic_sensitivity_models_ewas_chunks.R 500001 600000 ${i}
	Rscript Age22_epic_sensitivity_models_ewas_chunks.R 600001 700000 ${i}
	Rscript Age22_epic_sensitivity_models_ewas_chunks.R 700001 800000 ${i}
	Rscript Age22_epic_sensitivity_models_ewas_chunks.R 800001 900000 ${i}
done

cd <path to Age 22 sensitivity output>

for i in {1..7} {10..13};
do
	cat y22_epic_start_cpg_100000_model${i}.ewasoutput.txt > y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_200000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_300000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_400000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_500000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_600000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_700000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_800000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
awk '{if(NR>1){print $0}}' y22_epic_start_cpg_900000_model${i}.ewasoutput.txt >> y22_epic_model${i}.ewasoutput.txt
done
