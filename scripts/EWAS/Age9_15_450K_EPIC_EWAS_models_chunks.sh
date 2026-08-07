#!/bin/bash -l
# Author: Lisa Schneper
# Last update: Aug 6, 2026

# NEED TO SET FOLLOWING PATHS

SCRIPT_DIR="PATH_TO_SCRIPT_DIRECTORY"
Age9_OUTPUTDIR="PATH_TO_Age9_OUTPUT_DIRECTORY"
Age15_OUTPUTDIR="PATH_TO_Age15_OUTPUT_DIRECTORY"

##### RUN EWAS IN CHUNKS (faster on our machine)
module load R/4.5.1
cd ${SCRIPT_DIR}

###### Age 9
for i in 1 2 3 4 6 10 11; do
	Rscript Age9_450K_sensitivity_models_ewas_chunks_070226.R 1 100000 ${i}
	Rscript Age9_450K_sensitivity_models_ewas_chunks_070226.R 100001 200000 ${i}
	Rscript Age9_450K_sensitivity_models_ewas_chunks_070226.R 200001 300000 ${i}
	Rscript Age9_450K_sensitivity_models_ewas_chunks_070226.R 300001 400000 ${i}
	Rscript Age9_450K_sensitivity_models_ewas_chunks_070226.R 400001 500000 ${i}
	Rscript Age9_epic_sensitivity_models_ewas_chunks_070226.R 1 100000 ${i}
	Rscript Age9_epic_sensitivity_models_ewas_chunks_070226.R 100001 200000 ${i}
	Rscript Age9_epic_sensitivity_models_ewas_chunks_070226.R 200001 300000 ${i}
	Rscript Age9_epic_sensitivity_models_ewas_chunks_070226.R 300001 400000 ${i}
	Rscript Age9_epic_sensitivity_models_ewas_chunks_070226.R 400001 500000 ${i}
	Rscript Age9_epic_sensitivity_models_ewas_chunks_070226.R 500001 600000 ${i}
	Rscript Age9_epic_sensitivity_models_ewas_chunks_070226.R 600001 700000 ${i}
	Rscript Age9_epic_sensitivity_models_ewas_chunks_070226.R 700001 800000 ${i}
	Rscript Age9_epic_sensitivity_models_ewas_chunks_070226.R 800001 900000 ${i}
done

####### Age 15
for i in 1 2 3 4 5 6 7 10 11; do
	Rscript Age15_450K_sensitivity_models_ewas_chunks_070226.R 1 100000 ${i}
	Rscript Age15_450K_sensitivity_models_ewas_chunks_070226.R 100001 200000 ${i}
	Rscript Age15_450K_sensitivity_models_ewas_chunks_070226.R 200001 300000 ${i}
	Rscript Age15_450K_sensitivity_models_ewas_chunks_070226.R 300001 400000 ${i}
	Rscript Age15_450K_sensitivity_models_ewas_chunks_070226.R 400001 500000 ${i}
	Rscript Age15_epic_sensitivity_models_ewas_chunks_070226.R 1 100000 ${i}
	Rscript Age15_epic_sensitivity_models_ewas_chunks_070226.R 100001 200000 ${i}
	Rscript Age15_epic_sensitivity_models_ewas_chunks_070226.R 200001 300000 ${i}
	Rscript Age15_epic_sensitivity_models_ewas_chunks_070226.R 300001 400000 ${i}
	Rscript Age15_epic_sensitivity_models_ewas_chunks_070226.R 400001 500000 ${i}
	Rscript Age15_epic_sensitivity_models_ewas_chunks_070226.R 500001 600000 ${i}
	Rscript Age15_epic_sensitivity_models_ewas_chunks_070226.R 600001 700000 ${i}
	Rscript Age15_epic_sensitivity_models_ewas_chunks_070226.R 700001 800000 ${i}
	Rscript Age15_epic_sensitivity_models_ewas_chunks_070226.R 800001 900000 ${i}
done

# Age 9
cd ${Age9_OUTPUTDIR}
for model_num in 1 2 3 4 6 10 11; do
	cat y9_450K_sensitivity_1_100000_model${model_num}_model${model_num}.ewasoutput.txt >y9_450K_sensitivity_model${model_num}.ewasoutput.txt
	for i in {1..4}; do
		let "j=i+1"
		awk '{if(NR>1){print $0}}' y9_450K_sensitivity_${i}00001_${j}00000_model${model_num}_model${model_num}.ewasoutput.txt >>y9_450K_sensitivity_model${model_num}.ewasoutput.txt
	done
done

for model_num in 1 2 3 4 6 10 11; do
	cat y9_epic_sensitivity_1_100000_model${model_num}_model${model_num}.ewasoutput.txt >y9_epic_sensitivity_model${model_num}.ewasoutput.txt
	for i in {1..8}; do
		let "j=i+1"
		awk '{if(NR>1){print $0}}' y9_epic_sensitivity_${i}00001_${j}00000_model${model_num}_model${model_num}.ewasoutput.txt >>y9_epic_sensitivity_model${model_num}.ewasoutput.txt
	done
done

# Age 15
cd ${Age15_OUTPUTDIR}
for model_num in 1 2 3 4 5 6 7 10 11; do
	cat y15_450K_sensitivity_1_100000_model${model_num}_model${model_num}.ewasoutput.txt >y15_450K_sensitivity_model${model_num}.ewasoutput.txt
	for i in {1..4}; do
		let "j=i+1"
		awk '{if(NR>1){print $0}}' y15_450K_sensitivity_${i}00001_${j}00000_model${model_num}_model${model_num}.ewasoutput.txt >>y15_450K_sensitivity_model${model_num}.ewasoutput.txt
	done
done

for model_num in 1 2 3 4 5 6 7 10 11; do
	cat y15_epic_sensitivity_1_100000_model${model_num}_model${model_num}.ewasoutput.txt >y15_epic_sensitivity_model${model_num}.ewasoutput.txt
	for i in {1..8}; do
		let "j=i+1"
		awk '{if(NR>1){print $0}}' y15_epic_sensitivity_${i}00001_${j}00000_model${model_num}_model${model_num}.ewasoutput.txt >>y15_epic_sensitivity_model${model_num}.ewasoutput.txt
	done
done
