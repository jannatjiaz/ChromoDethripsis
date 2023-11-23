#!/bin/bash
#BSUB -J BAM_basecounter[1-10]
#BSUB -q normal
#BSUB -G team154-grp
#BSUB -e /lustre/scratch117/casm/team154/pc8/MetaNorm/Variant_calling/logs/basecounter_stderr_%I
#BSUB -o /lustre/scratch117/casm/team154/pc8/MetaNorm/Variant_calling/logs/basecounter_stdout_%I
#BSUB -n 1
#BSUB -M 1000
#BSUB -R "select[mem>1000] rusage[mem=1000] span[hosts=1]"

module load bgzip
module load tabix
log_dir=/lustre/scratch117/casm/team154/pc8/Jannat_pacbio/get_SNP_counts/logs/
output_dir=/lustre/scratch117/casm/team154/pc8/Jannat_pacbio/get_SNP_counts/Farm_output/

# Extract current sample
#LSB_JOBINDEX=1
echo $LSB_JOBINDEX
readarray positions < /lustre/scratch117/casm/team154/pc8/Jannat_pacbio/get_SNP_counts/Female_10Mb_bins_catchup.txt
pos="${positions[${LSB_JOBINDEX}-1]}"
pos=$(printf $pos)
echo $pos

sample_bam=/lustre/scratch117/casm/team154/pc8/Jannat_pacbio/get_SNP_counts/WTSI-OESO_152_2.sample.dupmarked.bam
echo $sample_bam
vcf_dir=/lustre/scratch117/casm/team154/ji2/OesophagealLongRead_WTSI-OESO_152/phaseCCS/whatshap/phased_variants/

julia /lustre/scratch117/casm/team154/pc8/Jannat_pacbio/get_SNP_counts/DRIVER_get_SNPs.jl $sample_bam $pos $vcf_dir $output_dir
