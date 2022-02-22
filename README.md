# CTResolver
## Haplotype resolution of chromothriptic and rearranged chromosomes 

Here we present a methodology for reconstructing chromothriptic and complex rearrangments which are restricted to a single allele. It will not work if the complex rearrangments affect both alleles but it is possible to reconstruct more complex rearrangments but this requires modification of the assignment rules. 

This analysis is done on a chromosome basis in order to make it easy to parallelizable and memory and time efficeint so split all bam and vcf files by chromosome. This is also important as haplotype resolved genome assemblies are done on a per chromosome basis.

## Whatshap phasing 

Whatshap can be found here: https://whatshap.readthedocs.io/en/latest/ 

In order run this code, you will need first need to run whatshap to aquire the phasing information. Although multiple sequencing modalitiies can be input into whatshap, we recommend that you only use CCS reads as this will lead to a high number of phase blocks but a low switching error rate. This is important as the presence of complex rearagements an subclonal regions may result in more switch errors than when phasing more normal samples. We recommend that the vcf used for phasing is filtered 

You will need to run: 

1) whatshap phase - to phase the variants into haplotype blocks
2) whatshap haplotag - to annotate the origianl CCS reads based on the SNPs that are present (input vcfs must be zipped and indexed)
3) whatshap stats - to the positions of the pahse blocks. 

You can split the haplotagged bam file into haplotype 1 and haplotyp2 using the HP ta, HP:i:1 and HP:i:2 respectively. For example for haplotype 1:
samtools view <haplotagged file> | awk -v tag="HP:i:1" '($0 ~ /^@/ || index($0,tag)>0)' | samtools view -b > <haplotype 1 bam> 
  

## Required input
  
This code is writing in python. The following modules need to be installed:

1) pysam 
2) pandas as pd  
  
### Resolving instances of chromothripsis with < 100 breakpoints 

If the denstiy of breakpoints on the chromothriptic region is not higher than the average density across the other samples, ou will not be able to confidently assign a read to the chromothriptic over the wild-type. Therefore SV phasing information will not be used. 

In order run this code, :

1) haplotagged bam file - output from whatshap haplotag
2) txt file of all read names for haplotype 1, haplotype 2 and the unresolved reads
3) a tab seperated file of haplotype 2 reads with columns <readname> <chr> <start_position>
4) a tab seperated file of haplotype 2 reads with columns <readname> <chr> <start_position>
5) whathap block information - output from whatshap stat
6) a phased haplotype file with colmns <chr> <Haplotype start> <Number of SNPs> <ave_VAF_block1> <ave_VAF_block2> <ave_depth_block1> <ave_depth_block1>.
7) region of loh file with columns <chr> <start> <end>

With these files you can run CTResolveLowDensity.py:
  

python phase_ccs_read.py --chromosome ${CHR} \
  --chromosome <chromosome>
  --bam <path to haplotagged bamfile (1) >
  --readnames <path to readnames list (2) > \
  --haplotype1_reads <readnames pos chr haplotype1 (3) > \
  --haplotype2_reads <readnames pos chr haplotype1 (4) > \
  --haplotype_blocks  <path to the haplotype blocks produced by whatshap stats (5) > \
  --phased_haplotype_blocks <path to the phased haplotype blocks produced by Peters code (6)>  \
  --loh_regions <path to loh file (7) >
  --outputdir <path to the output directory>

This will produce sam files that contain reads that phase to haplotype 1 (called haplotype1_reads_all_<chromosome>_ccs.sam) and haplotype 1 (called haplotype1_reads_all_<chromosome>_ccs.sam).
  
  
### Resolving instances of chromothripsis with > 100 breakpoints 

If the denstiy of breakpoints on the chromothriptic region is higher than the average density across the other samples, it will be informative to use SV presence to assign haplotype blocks to one haplotype over the other. 

In order to do this run a long read SV caller which will output readnames with the SV calls 



