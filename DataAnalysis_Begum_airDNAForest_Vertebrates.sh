#######################################
#    Data Analyses for manuscript: Airborne environmental DNA captures terrestrial vertebrate diversity in nature
# date: 25/04/22
# Analysis of data from 12Svert (Riaz)
# Christina Lynggaard
#######################################

#To obtain raw sequence data go to: https://erda.ku.dk/archives/edc5a1a0a89ff2f9ece0ba7ee12f3196/published-archive.html

##################  1. MAKE ALL THE FILES READY ####################

#rename the files so they are easier to work with:
mv 12Sv_PSinfo_DNAir_Aa_Field_replicate1.txt 12Sv_PSinfo_DNAir_Aa.txt #contains tag information for field replicate 1 samples
mv 12Sv_PSinfo_DNAir_Aa_Field_replicate2.txt 12SvR2_PSinfo_DNAir_Aa.txt #contains tag information for field replicate 1 samples

# In the working folder the following 4 files should exist:
			12Sv_poolInfo.txt #name of pool(Column1) path to R1 fastq file(Column2) path to R2 fastw file(Colum 3) if not merged / can be .gz 
								#name of pool(Column1) path to the collapsed and merged files 
				#pool file is created after adapter removal is used
			12SvR2_poolInfo.txt #For Field Replicate 2 samples.
			12Sv_Primers.txt # just 2 columns: sequence of primerF and primerR 
			12Sv_PSinfo_DNAir_Aa.txt  #name of sample(Column1) Tag#(Fwd/Column2) Tag#(Rvrs/Column3) Pool#(Column4; need to be same name as in poolInfo.txt)
			12SvR2_PSinfo_DNAir_Aa.txt
			12Sv_Tags.txt # tag name(Column1) tag sequence(Column2)

##################  2. Runc FASTQC ################## 
module load fastqc/v0.11.9

#On the first field replicate samples
fastqc Riaz1Aamosen_S1_L001_R1_001.fastq.gz 
fastqc Riaz1Aamosen_S1_L001_R2_001.fastq.gz 

fastqc Riaz2Aamosen_S2_L001_R1_001.fastq.gz
fastqc Riaz2Aamosen_S2_L001_R2_001.fastq.gz

fastqc Riaz3Aamosen_S3_L001_R1_001.fastq.gz
fastqc Riaz3Aamosen_S3_L001_R2_001.fastq.gz

fastqc Riaz4Aamosen_S4_L001_R1_001.fastq.gz
fastqc Riaz4Aamosen_S4_L001_R2_001.fastq.gz

#Now on the second field replicate samples
fastqc Riaz1AamosenR2_S1_L001_R1_001.fastq.gz
fastqc Riaz1AamosenR2_S1_L001_R2_001.fastq.gz

fastqc Riaz2AamosenR2_S2_L001_R1_001.fastq.gz
fastqc Riaz2AamosenR2_S2_L001_R2_001.fastq.gz

fastqc Riaz3AamosenR2_S3_L001_R1_001.fastq.gz
fastqc Riaz3AamosenR2_S3_L001_R2_001.fastq.gz

fastqc Riaz4AamosenR2_S4_L001_R1_001.fastq.gz
fastqc Riaz4AamosenR2_S4_L001_R2_001.fastq.gz

##### 3. ADAPTOR REMOVAL 

#Load AdaptorRemoval
module load AdapterRemoval/v2.3.1  

AdapterRemoval --file1 Riaz1Aamosen_S1_L001_R1_001.fastq.gz --file2 Riaz1Aamosen_S1_L001_R2_001.fastq.gz --minlength 100 --maxlength 160 --shift 5 --basename pool1_merged --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse
AdapterRemoval --file1 Riaz2Aamosen_S2_L001_R1_001.fastq.gz --file2 Riaz2Aamosen_S2_L001_R2_001.fastq.gz --minlength 100 --maxlength 160 --shift 5 --basename pool2_merged --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse
AdapterRemoval --file1 Riaz3Aamosen_S3_L001_R1_001.fastq.gz --file2 Riaz3Aamosen_S3_L001_R2_001.fastq.gz --minlength 100 --maxlength 160 --shift 5 --basename pool3_merged --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse
AdapterRemoval --file1 Riaz4Aamosen_S4_L001_R1_001.fastq.gz --file2 Riaz4Aamosen_S4_L001_R2_001.fastq.gz --minlength 100 --maxlength 160 --shift 5 --basename pool4_merged --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse
AdapterRemoval --file1 Riaz1AamosenR2_S1_L001_R1_001.fastq.gz --file2 Riaz1AamosenR2_S1_L001_R2_001.fastq.gz --minlength 100 --maxlength 160 --shift 5 --basename pool1R2_merged --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse
AdapterRemoval --file1 Riaz2AamosenR2_S2_L001_R1_001.fastq.gz --file2 Riaz2AamosenR2_S2_L001_R2_001.fastq.gz --minlength 100 --maxlength 160 --shift 5 --basename pool2R2_merged --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse
AdapterRemoval --file1 Riaz3AamosenR2_S3_L001_R1_001.fastq.gz --file2 Riaz3AamosenR2_S3_L001_R2_001.fastq.gz --minlength 100 --maxlength 160 --shift 5 --basename pool3R2_merged --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse
AdapterRemoval --file1 Riaz4AamosenR2_S4_L001_R1_001.fastq.gz --file2 Riaz4AamosenR2_S4_L001_R2_001.fastq.gz --minlength 100 --maxlength 160 --shift 5 --basename ppool4R2_merged --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse


#	Concatenate the merged collapsed files 
cat pool1_merged.collapsed pool1_merged.collapsed.truncated > pool1_merged.fastq
cat pool2_merged.collapsed pool2_merged.collapsed.truncated > pool2_merged.fastq
cat pool3_merged.collapsed pool3_merged.collapsed.truncated > pool3_merged.fastq
cat pool4_merged.collapsed pool4_merged.collapsed.truncated > pool4_merged.fastq
#OBS! the 12Sv_poolInfo.txt file will have the path to where these pool*.fastq files are

cat pool1_merged.collapsed pool1R2_merged.collapsed.truncated > pool1R2_merged.fastq
cat pool2_merged.collapsed pool2R2_merged.collapsed.truncated > pool2R2_merged.fastq
cat pool3_merged.collapsed pool3R2_merged.collapsed.truncated > pool3R2_merged.fastq
cat pool4_merged.collapsed pool4R2_merged.collapsed.truncated > pool4R2_merged.fastq
#OBS! the 12SvR2_poolInfo.txt file will have the path to where these pool*R2.fastq files are


#check the quality after merging
fastqc pool1_merged.fastq
fastqc pool2_merged.fastq 
fastqc pool3_merged.fastq
fastqc pool4_merged.fastq 
fastqc pool1R2_merged.fastq
fastqc pool2R2_merged.fastq 
fastqc pool3R2_merged.fastq
fastqc pool4R2_merged.fastq 


################## Begum ############################
##### 4. 'sort' Files with Begum  ###                                         
		
module load python/v2.7.12 
module load Begum

mkdir begum_2mismatches
Begum sort -p 12Sv_Primers.txt -t 12Sv_Tags.txt -s 12Sv_PSinfo_DNAir_Aa.txt -l 12Sv_poolInfo.txt -pm 2 -d begum_2mismatches -o begum_12Sv
Begum sort -p 12Sv_Primers.txt -t 12Sv_Tags.txt -s 12SvR2_PSinfo_DNAir_Aa.txt -l 12SvR2_poolInfo.txt -pm 2 -d begum_2mismatches -o begum_12SvR2

cd begum_2_mismatches

##### 5. 'filter' Files with Begum  ### 
# Select sequences found in 1 out of 4 PCR replicates

mkdir filtered_1out4
Begum filter -i begum_12Sv -s 12Sv_PSinfo_DNAir_Aa.txt -p 0.25 -m 35 -l 90 -d filtered_1out4 -o filtered_025_35_90
Begum filter -i begum_12Sv -s 12SvR2_PSinfo_DNAir_Aa.txt -p 0.25 -m 35 -l 90 -d filtered_1out4 -o filteredR2_025_35_90

#explore the negative and positive controls
grep bl filtered_025_35_90.fna 
grep POS filtered_025_35_90.fna
cd ..

#Combine filtered data from both field replicates into one file
cat filtered_025_35_90.fna filteredR2_025_35_90 > combined_filtered_025_35_90.fna


################## OTU CLUSTERING ############################
module load sumaclust

#Convert working-file to required format.
#Get the python scripts from here: https://github.com/shyamsg/DAMe/tree/master/bin
python convertToUSearch.py -i combined_filtered_025_35_90.fna -lmin 90 -lmax 110

# Create OTU clusters with sumaclust 
sumaclust -e FilteredReads.forsumaclust.fna -t 0.97 -F combined_DNAirAa_Riaz_OTUs_sumaclust.fna

#Create working files
python tabulateSumaclust.py -i combined_DNAirAa_Riaz_OTUs_sumaclust.fna -o DNAirAa_Combined_Riaz_025_min35.txt -blast

#Download files for further processing