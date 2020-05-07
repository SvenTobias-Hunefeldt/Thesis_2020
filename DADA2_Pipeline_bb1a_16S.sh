#Make output directories path
WorkingPATH="/SCRATCH/Sven_Data/DADA2/Biofilm/16S"

#Need to pre-make samples folder and insert metadata file
mkdir ${WorkingPATH}/test_sample_IDs
mkdir ${WorkingPATH}/test_sample_IDs/R1/
mkdir ${WorkingPATH}/test_sample_IDs/R2/

mkdir ${WorkingPATH}/Samples

mv ${WorkingPATH}/taxa_DB ${WorkingPATH}/Samples/ #Maybe copy? Depends on if you run into issues I guess.


cp ${WorkingPATH}/*.txt ${WorkingPATH}/Samples/Mapping_file.txt


#Demultiplex R1
split_libraries_fastq.py -i ${WorkingPATH}/Undetermined_S0_L001_R1_001.fastq -b ${WorkingPATH}/Undetermined_S0_L001_I1_001.fastq -m ${WorkingPATH}/*.txt -n 0 -o ${WorkingPATH}/test_demultplx_seqs_R1.fna --rev_comp_mapping_barcodes --rev_comp_barcode --store_demultiplexed_fastq &

#Demultiplex R2
split_libraries_fastq.py -i ${WorkingPATH}/Undetermined_S0_L001_R2_001.fastq -b ${WorkingPATH}/Undetermined_S0_L001_I1_001.fastq -m ${WorkingPATH}/*.txt -n 0 -o ${WorkingPATH}/test_demultplx_seqs_R2.fna --rev_comp_mapping_barcodes --rev_comp_barcode --store_demultiplexed_fastq

wait $(jobs -p)


#Split R1 by sample ID
split_sequence_file_on_sample_ids.py -i ${WorkingPATH}/test_demultplx_seqs_R1.fna/seqs.fastq --file_type fastq -o ${WorkingPATH}/test_sample_IDs/R1/ &


#Split R2 by sample ID
split_sequence_file_on_sample_ids.py -i ${WorkingPATH}/test_demultplx_seqs_R2.fna/seqs.fastq --file_type fastq -o ${WorkingPATH}/test_sample_IDs/R2/

wait $(jobs -p)




#Change path for file renaming
cd ${WorkingPATH}/test_sample_IDs/R1/

#Rename files
for file in *.fastq
 do
    mv "${file}" "${WorkingPATH}/Samples/${file}_R1_001.fastq"
 done 
 
#Change path for file renaming
 cd ${WorkingPATH}/test_sample_IDs/R2/

#Rename files
for file in *.fastq
 do
    mv "${file}" "${WorkingPATH}/Samples/${file}_R2_001.fastq"
 done 



#Run R script
Rscript ${WorkingPATH}/DADA2_16S.R








