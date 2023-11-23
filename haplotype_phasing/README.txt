bsub -J BAM_basecounter[1-294] -q normal -G team154-grp -e /lustre/scratch117/casm/team154/pc8/Jannat_pacbio/get_SNP_counts/logs/basecounter_stderr_%I -o /lustre/scratch117/casm/team154/pc8/Jannat_pacbio/get_SNP_counts/logs/basecounter_stdout_%I -n 1 -M 1000 -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" /lustre/scratch117/casm/team154/pc8/Jannat_pacbio/get_SNP_counts/Run_basecounter.sh

# Then, in Julia:

#!/usr/bin/julia

include("Haplotype_functions.jl")
ENV["GKSwstype"]="100"

sampleID = "OESO_117"

chromosomes = "chr" .* [string.(1:22);]
snp_dir = "/lustre/scratch117/casm/team154/pc8/Jannat_pacbio/get_SNP_counts/$sampleID/"
snp_files = readdir(snp_dir)
for i in chromosomes
    println(i)
    sub_files = snp_dir .* snp_files[occursin.(Regex("_$(i)_"), snp_files)]
    write("$(i)_hetSNP_cts.txt", read(pipeline(`awk 'FNR==1 && NR!=1{next;}{print}' $sub_files`, 
                `sort -n +1`)))
    
    # Call the phasing algorithm to generate DataFrame per SNP of haplotype information
    a = phasing_het_SNPs("/lustre/scratch117/casm/team154/ji2/OesophagealLongRead_WTSI-$sampleID/phaseCCS/whatshap/phased_variants/variants_phased_nosimplerepeats_filtered_ccs_$i.vcf", 
                Interval(i, 1, 260000000), "$(i)_hetSNP_cts.txt", "$(i)_full_chr.pdf")

    # Summarise the information per haplotype block
    b = allele_fraction_hap_blocks(a)
    b = hcat(DataFrame(Chr = repeat([i;], size(b,1))), b)

    # Write to file
    CSV.write("Haplotype_blocks_$i.txt", b, delim="\t")
end

hap_files = "Haplotype_blocks_" .* chromosomes .* ".txt"
ct_files = chromosomes .* "_hetSNP_cts.txt"
write("$(sampleID)_haplotype_blocks.txt", read(pipeline(`awk 'FNR==1 && NR!=1{next;}{print}' $hap_files`)))
rm.(hap_files)
rm.(ct_files)


