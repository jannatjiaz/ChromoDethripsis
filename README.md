# CTResolver
## Haplotype resolution of chromothriptic and rearranged chromosomes 

Here we present a methodology for reconstructing chromothriptic and complex rearrangments which are restricted to a single allele. It will not work if the complex rearrangments affect both alleles but it is possible to reconstruct more complex rearrangments but this requires modification of the assignment rules. 

This analysis is doen on a chromosome basis in order to make it easy to parallelizable and memory and time efficeint so split all bam and vcf files by chromosome.

## Whatshap phasing 

Whatshap can be found here: https://whatshap.readthedocs.io/en/latest/ 

In order run this code, you will need first need to run whatshap to aquire the phasing information. Although multiple sequencing modalitiies can be input into whatshap, we recommend that you only use CCS reads as this will lead to a high number of phase blocks but a low switching error rate. This is important as the presence of complex rearagements an subclonal regions may result in more switch errors than when phasing more normal samples. 

You will need to run: 

1) whatshap phase - to phase the variants into haplotype blocks
2) whatshap haplotag - to annotate the origianl CCS reads based on the SNPs that are present (input vcfs must be zipped and indexed)
3) whatshap stats - to the positions of the pahse blocks. 

You can split the haplotagged bam file into haplotype 1 and haplotyp2 using the HP ta, HP:i:1 and HP:i:2 respectively. For example for haplotype 1:
samtools view <haplotagged file> | awk -v tag="HP:i:1" '($0 ~ /^@/ || index($0,tag)>0)' | samtools view -b > <haplotype 1 bam> 
  

## Required input

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

With these files you can run CTResolveLowDensity:
  


### Resolving instances of chromothripsis with > 100 breakpoints 

If the denstiy of breakpoints on the chromothriptic region is higher than the average density across the other samples, it will be informative to use SV presence to 



