# ChromoDethripsis
## Haplotype resolution of chromothriptic and rearranged chromosomes 

Here we present a methodology for reconstructing chromothriptic and complex rearrangments which are restricted to a single allele. It will not work if the complex rearrangments affect both alleles but it is possible to reconstruct more complex rearrangments but this requires modification of the assignment rules. 

This analysis is done on a chromosome basis in order to make it easy to parallelizable and memory and time efficeint so split all bam and vcf files by chromosome. This is also important as haplotype resolved genome assemblies are done on a per chromosome basis.

## Whatshap phasing 

Whatshap can be found here: https://whatshap.readthedocs.io/en/latest/ 

In order run this code, you will need first need to run whatshap to aquire the phasing information. Although multiple sequencing modalitiies can be input into whatshap, we recommend that you only use CCS reads as this will lead to a high number of phase blocks but a low switching error rate. This is important as the presence of complex rearagements an subclonal regions may result in more switch errors than when phasing more normal samples. We recommend that the vcf used for phasing is filtered to remove repeat regions or low complexity regions as annotated by UCSC. This is to refuce switch errors.

You will need to run: 

1) whatshap phase - to phase the variants into haplotype blocks
2) whatshap haplotag - to annotate the origianl CCS reads based on the SNPs that are present (input vcfs must be zipped and indexed)
3) whatshap stats - to the positions of the pahse blocks. 

You can split the haplotagged bam file into haplotype 1 and haplotype 2 using the HP tag, HP:i:1 and HP:i:2 respectively. For example for haplotype 1:
> samtools view input_haplotagged.bam | awk -v tag="HP:i:1" '($0 ~ /^@/ || index($0,tag)>0)' | samtools view -b > output_haplotype1.bam
  

## CCS reads haplotype resolution
  
This code is writing in python. You must use python 3. The following modules need to be installed:

1) pysam 
2) pandas
3) xlrd==1.2.0
  
### Resolving instances of chromothripsis with < 100 breakpoints 

If the denstiy of breakpoints on the chromothriptic region is not higher than the average density across the other samples, ou will not be able to confidently assign a read to the chromothriptic over the wild-type. Therefore SV phasing information will not be used. 

Files needed in order run this code are:

1) haplotagged bam file - output from whatshap haplotag
2) txt file of all read names for haplotype 1, haplotype 2 and the unresolved reads
3) a tab seperated file of haplotype 2 reads with 3 columns: readname, chr, start_position
4) a tab seperated file of haplotype 2 reads with 3 columns: readname, chr, start_position
5) whathap block information - output from whatshap stat
6) a phased haplotype file with 7 columns: chr, haplotype, start, number of SNPs, ave_VAF_block1, ave_VAF_block2, ave_depth_block1, ave_depth_block1
7) region of loh file with 3 columns: chr start end

With these files you can run CTResolveLowDensity.py. If the copy number of the chromothriptic and wild-type chromosomes are equal (i.e. 1:1) or if the copy number of the chromothriptic chromosome is higher than the copy number of the wild-type chromsomse, use CTResolveLowDensity.py/ If the copy number of the wild-type chromosome is higher than the chromothriptic chromsomse, use CTResolveLowDensityDupWildtype.py. The parameters for both scripts are the same, only the assignment rules change. Parameters include: 

> python CTResolveLowDensity.py \
>  --chromosome <chromosome> \
>  --bam <path to haplotagged bamfile (1) > \
>  --readnames <path to readnames list (2) > \
>  --haplotype1_reads <readnames pos chr haplotype1 (3) > \
>  --haplotype2_reads <readnames pos chr haplotype1 (4) > \
>  --haplotype_blocks  <path to the haplotype blocks produced by whatshap stats (5) > \
>  --phased_haplotype_blocks <path to the phased haplotype blocks produced by Peters code (6)>  \
>  --loh_regions <path to loh file (7) > \
>  --outputdir <path to the output directory> \

This will produce sam files that contain reads that phase to haplotype 1 (called haplotype1_reads_all_<chromosome>_ccs.sam) and haplotype 1 (called haplotype1_reads_all_<chromosome>_ccs.sam).
  
  
### Resolving instances of chromothripsis with > 100 breakpoints 

If the denstiy of breakpoints on the chromothriptic region is higher than the average density across the other samples, it will be informative to use SV presence to assign haplotype blocks to one haplotype over the other. At present this is only suitable for human chromosomes mapped to GRCH38. 

In order to do this, run a long read SV caller which will output readnames with the SV calls. Examples include Sniffles https://github.com/fritzsedlazeck/Sniffles or pbsv https://github.com/PacificBiosciences/pbsv but any caller can be used. In order to make the variant caller bedpe compatible with the resolution code, ensure that columns are: CHR1, POS1, POS2, CHR2, POS3, POS4, ID, QUAL, STRAND1, STRAND2, REF, TYPE, FILTER, INFO, GT:DR:DV, INFO2. This is the default output from Sniffles. 

Files needed in order run this code are:

1) haplotagged bam file - output from whatshap haplotag
2) bedpe file - structural variants called across all chromosomes not just chromosome of interest (used to determine baseline SV rate)
3) heterozygous snps file - vcf file used to for phasing by whatshap. The column format should match strelka2 and header should be removed
4) txt file of all read names for haplotype 1, haplotype 2 and the unresolved reads
5) a tab seperated file of haplotype 2 reads with 3 columns: readname, chr, start_position
6) a tab seperated file of haplotype 2 reads with 3 columns: readname, chr, start_position
7) whathap block information - output from whatshap stat
8) phasing haplotype blocks file - information about features of blocks produced by peters code with 7 columns: Haplotype, Number_SNPs, ave_VAF_block1, ave_depth_block1, ave_depth_block2
9) phased VCF file - phased variants produced by whatshap
10) region of loh file with 3 columns: chr start end
  

With these files you can run CTResolveLowDensity.py. It assumes that there are equal copies of the chromothrpitic and wild-type chromosome. Parameters include: 

> python CTResolveHighDensity.py \
>  --chromosome <chromosome> \
>  --bam <path to haplotagged bamfile (1) > \
>  --bedpe <path to SV bedpe (2) > \
>  --heterozygousSnps <path to snps used for phasing (3)> \
>  --readnames <path to readnames list (4) > \
>  --haplotype1_reads <readnames pos chr haplotype1 (5) > \
>  --haplotype2_reads <readnames pos chr haplotype1 (6) > \
>  --haplotype_blocks  <path to the haplotype blocks produced by whatshap stats (7) > \
>  --phased_haplotype_blocks <path to the phased haplotype blocks produced by Peters code (8)>  \
>  --phased_vcf <phased variants produced by whatshap (9) >   \
>  --loh_regions <path to loh file (10) > \
>  --outputdir <path to the output directory> \

### Next steps
  
ChomoDethripsis will produce a bam file for each haplotype. From this conventional assembly methods can be ultilised. It is important to remember that resolution of complex strucutral variation will require manual checking of results to ensure main features of chromothrpisis have been accurately reconstructed. 
  
  
  

  
  



