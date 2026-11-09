Author: Eve Schoeffler
Date: July 26, 2024

For small datasets, use whole_ENmix_preprocessing.R

For large datasets, use the scripts as outlined below.

######################    

1. Split large RGChannelSet into plate-sized batches
QCinfo should be run on an RGChannelSet to find bad samples and probes. If the RGSet is large (over 2000 samples), QCinfo can't be run directly. Instead. the RGSet must be split into plate-sized batches (~96, though this isn't a strict number). Decide the number of batches that you want and write it in the script. Use this script:
  
 split_combinedRGSet_intoBatches_forQCsplitter.R
 
######################   

2. Run QCinfo on batches
After the plate-sized batches have saved, you will run a version of QCinfo on each of the batches. For the number of batches you created, run the following script:

  for batch in $(seq 1 num_batches); do Rscript QCsplitter_perBatch.R ${batch}; done
    (num_batches needs to be set by you)
   
######################     

3. Combine QC chunks
Make a folder called Final_QC in the same folder as your large RGSet. Combine each of the QC batches with the following script. It will save a final QCinfo object for the whole dataset. You need to go in, change the paths, and the number of batches (plates).

  combine_EPIC_QCinfo_final.R

######################      

4. Split RGSet into batches of ~800 samples
Split your RGSet into smaller batches that are saved in Preprocess/RGSets. Rather than ~96 samples, now make it into ~800 samples. Change paths and edit the following script:

  split_bigRGSet_forPreprocessing.R
  
######################    

5. Preprocessing
Run the preprocessing script on each batch! Edit the following to correct the paths to your needs. Run one section at a time (comment out other sections) to keep memory use low. If you're running on sections of ~800 samples, preprocessENmix (the first step) takes about 12 hours. The next steps are quicker.

  for batch in $(seq 1 num_batches); do Rscript qc_preprocess_batch.R ${batch}; done
  
######################    

6. Normalization & Dye Bias Correction
Combine the ~800 sample methyl sets created in the previous step. Then run quantile normalization and rcp to produce the final beta matrix.

  preprocess_last_steps.R

  


