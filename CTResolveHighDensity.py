#! /usr/bin/python

import pysam 
import pandas as pd
import random
import os 
from os import system
import argparse 

parser = argparse.ArgumentParser(description='usage:  CTResolveHighDensity.py --chromosome chromosome --bam bam --bedpe bedpe --hetsnps hetsnps --readnames readnames\
	--haplotype1_reads haplotype1_reads --haplotype2_reads haplotype2_reads --haplotype_blocks haplotype_blocks --phased_haplotype_blocks phased_haplotype_blocks \
	--phased_vcf phased_vcf --loh_regions --loh_regions --outputdir path')
parser.add_argument('--chromosome', help='chromosome',required=True)
parser.add_argument('--bam', help='path to bamfile',required=True)
parser.add_argument('--bedpe', help='path to bedpe',required=True)
parser.add_argument('--heterozygousSnps', help='path to heterozygous snp file',required=True)
parser.add_argument('--readnames', help='path to readnames list',required=True)
parser.add_argument('--haplotype1_reads', help='readnames pos chr haplotype1 ccs',required=True)
parser.add_argument('--haplotype2_reads', help='readnamespos chr  haplotype2 ccs',required=True)
parser.add_argument('--haplotype_blocks', help='path to the haplotype blocks produced by whatshap stats',required=True)
parser.add_argument('--phased_haplotype_blocks', help='path to the phased haplotype blocks produced by Peters code',required=True)
parser.add_argument('--phased_vcf', help='path to phased vcf file from whatshap',required=True)
parser.add_argument('--loh_regions', help='path to loh regions file',required=True)
parser.add_argument('--outputdir', help='path to the output directory',required=True)
args = parser.parse_args()


chromosome=args.chromosome
bamfile=args.bam
bedpefile=args.bedpe
hetsnps=args.heterozygousSnps
readnames_file=args.readnames
haplotype1_readnames_chr_pos_file=args.haplotype1_reads
haplotype2_readnames_chr_pos_file=args.haplotype2_reads
haplotype_blocks_file=args.haplotype_blocks
phasing_haplotype_blocksfile=args.phased_haplotype_blocks
phased_vcf_file = args.phased_vcf
loh_regions_file=args.loh_regions
outputdir=args.outputdir


bamfile="/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/python_phasingSVs/haplotype_phasing/ccs_no_simple_repeats/haplotag_chr6_phased_ccs_no_simple_repeats.bam"
bedpefile =  "/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/python_phasingSVs/sniffles/OESO_103_sniffles_hg38.bedpe"
hetsnps = "/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/python_phasingSVs/snp_phasing/variants_chr6_allhet_no_simple_repeats_filtered_nohead.vcf"
readnames_file = "/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/python_phasingSVs/haplotype_phasing/ccs_no_simple_repeats/chr6_all_read_names.txt"
haplotype1_readnames_chr_pos_file = "/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/python_phasingSVs/haplotype_phasing/ccs_no_simple_repeats/haplotype1_read_chr_pos_new.txt"
haplotype2_readnames_chr_pos_file = "/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/python_phasingSVs/haplotype_phasing/ccs_no_simple_repeats/haplotype2_read_chr_pos_new.txt"
haplotype_blocks_file = "/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/python_phasingSVs/haplotype_phasing/ccs_no_simple_repeats/ccs_blocks.tsv"
phasing_haplotype_blocksfile = "Haplotype_blocks_chr6.txt"
phased_vcf_file = "/lustre/scratch117/casm/team154/ji2/OesophagealLongRead/python_phasingSVs/haplotype_phasing/ccs_no_simple_repeats/variants_phased_nosimplerepeats_filtered_ccs.vcf"
loh_regions = "tmp2"
outputdir = "phased_CCS_reads"


##LOAD DATA
#bam
bam = pysam.AlignmentFile(bamfile, "rb") #load bam file

#load SV calls from sniffles
bedpe = pd.read_csv(bedpefile, sep="\t", header=None ) #load bedpe and give it column names 
bedpe.columns = ["CHR1","POS1","POS2","CHR2","POS3","POS4","ID","QUAL","STRAND1","STRAND2","REF","TYPE","FILTER","INFO","GT:DR:DV","INFO2"]
SV_ID = bedpe["CHR1"].astype(str)+'_'+bedpe["POS2"].astype(str)+'_'+bedpe["STRAND1"]+'_'+bedpe["CHR2"].astype(str)+'_'+bedpe["POS4"].astype(str)+'_'+bedpe["STRAND2"]
SV_ID = SV_ID.tolist() #convert it to a list so you can use it as keys 
#make the SV dictionary where each key is an SV id and each key contains, SV_id, chr1, pos1, strand1, chr2, pos2, strand2 and reads
SV_dict={x:{} for x in SV_ID} #make an empty dictionary for each of the SVs
for item in range(len(SV_ID)):
	SV_dict[SV_ID[item]]= {"SV_id":SV_ID[item], "chrom1":bedpe["CHR1"][item], "pos1":bedpe["POS2"][item],"strand1":bedpe["STRAND1"][item],"chrom2":bedpe["CHR2"][item], "pos2":bedpe["POS4"][item],"strand2":bedpe["STRAND2"][item]}

reads_with_SVs_tmp = []
for SV in range(len(bedpe)):
	reads=bedpe["INFO"][SV].split(';')[10][7:].split(",")
	reads_with_SVs_tmp.append(reads)

reads_with_SVs = [y for x in reads_with_SVs_tmp for y in x]
reads_with_SVs= set(reads_with_SVs)

#load heterozygous snps 
snps = pd.read_csv(hetsnps, sep="\t", header=None )
snps.columns = ["CHR","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE"]
snp_ids = list(snps["CHR"].astype(str)+'_'+snps["POS"].astype(str))
snp_dict = {x:{} for x in snp_ids}
for snp in range(len(snp_ids)):
	snp_dict[snp_ids[snp]]={"snp_id":snp_ids[snp], "chr":snps["CHR"][snp], "pos":snps["POS"][snp], "ref":snps["REF"][snp], "alt":snps["ALT"][snp],"sample_info":snps["SAMPLE"][snp]}

#load readnames
readnames = pd.read_csv(readnames_file, sep="\t", header=None ) #load bedpe and give it column names 
readnames.columns=["reads"]
reads_dict={x:{} for x in list(readnames["reads"])} 

#haplotype1
haplotype1_read = pd.read_csv(haplotype1_readnames_chr_pos_file, sep="\t", header=None ) #load bedpe and give it column names 
haplotype1_read.columns = ['read', 'chr', 'pos']

#haplotype2
haplotype2_read = pd.read_csv(haplotype2_readnames_chr_pos_file, sep="\t", header=None ) #load bedpe and give it column names 
haplotype2_read.columns = ['read', 'chr', 'pos']

#haplotype blocks
haplotype_blocks = pd.read_csv(haplotype_blocks_file, sep="\t", skiprows=1, header=None) 
haplotype_blocks.columns = ['sample','chromosome','phase_set','from','to','variants']
#record block length
block_length =list(haplotype_blocks['to']-haplotype_blocks['from'])
#make a dictionary so the data is more accessible
haplotype_blocks_dict={x:{} for x in list(haplotype_blocks['phase_set'])} #make an empty dictionary for each of the SVs
for key in range(len(haplotype_blocks_dict)):
	haplotype_blocks_dict[haplotype_blocks['phase_set'][key]]={'chr':haplotype_blocks['chromosome'][key], 'from':haplotype_blocks['from'][key], 'to':haplotype_blocks['to'][key], 'length':block_length[key], 'variants':haplotype_blocks['variants'][key]}

#phasing
#phasing = pd.read_csv(phasing_haplotype_blocksfile, sep="\t", skiprows=1, header=None)
#phasing.columns = ["Haplotype","Number_SNPs","ave_VAF_block1","ave_depth_block1","ave_depth_block2"]
phasing = pd.read_excel(phasing_haplotype_blocksfile)

#make a dictionary so the data is more accessible
phasing_dict={x:{} for x in list(phasing['Haplotype'])} #make an empty dictionary for each of the SVs
for key in range(len(phasing_dict)):
	phasing_dict[phasing['Haplotype'][key]]={'Number_SNPs':phasing['Number_SNPs'][key], 'Number_LOH_Block_1':phasing['Number_LOH_Block_1'][key], 'Number_LOH_Block_2':phasing['Number_LOH_Block_2'][key], 'Number_SVs_Block_1':phasing['Number_SVs_Block_1'][key], 'Number_SVs_Block_2':phasing['Number_SVs_Block_2'][key]}

#find the block length, start and end and then match them with the phasing info
for key in phasing_dict:
	for haplotype in haplotype_blocks_dict:
		if key == haplotype:
			phasing_dict[key]['block_length']=haplotype_blocks_dict[int(key)]['length']
			phasing_dict[key]['from']=haplotype_blocks_dict[int(key)]['from']
			phasing_dict[key]['to']=haplotype_blocks_dict[int(key)]['to']

#load in the phased vcf:
phased_vcf = pd.read_csv(phased_vcf_file, sep='\t', header=None, skiprows=3)
phased_vcf.columns = ["CHR","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE"]
phased_vcf_ids = list(phased_vcf["CHR"].astype(str)+'_'+phased_vcf["POS"].astype(str))
phased_snp_dict = {x:{} for x in phased_vcf_ids}
for snp in range(len(phased_vcf_ids)):
	phased_snp_dict[phased_vcf_ids[snp]]={"snp_id":phased_vcf_ids[snp], "chr":phased_vcf["CHR"][snp], "pos":phased_vcf["POS"][snp], "ref":phased_vcf["REF"][snp], "alt":phased_vcf["ALT"][snp],"sample_info":phased_vcf["SAMPLE"][snp], "phasing":phased_vcf["SAMPLE"][snp].split(':')[0]}


#determine the average frequency of SVs on non-chromothriptic chromosomes - divide by 2 as there are two alleles
#maks the list of chromosomes 1 to 22
chr_list = []
for s in range(1, 23):
	chr_temp = 'chr'+str(s)
	chr_list.append(chr_temp)

chrs ={x:0 for x in chr_list} #make an empty dictionary for each of the chr
#cont the number of SVs on each chromosome
for sv in SV_dict:
	for chr in chrs:
		if SV_dict[sv]['chrom1']==chr:
			chrs[chr]=chrs[chr]+1

chomothriptic_chr = chrs[chromosome] # get the SV counts for the chromothriptic chromosome
del chrs[chromosome] #remove chr6

chr_lengths = [248946422,242183529,198285559,190204555,181528259,170795979,159335973,145128636,138384717,133787422,135076622,133265309,114354328,107033718,101981189,90328345,83247441,80363285,58607616,64434167,46699983,50808468] #for chr 1-22
chrs_len ={x:0 for x in chr_list} #make an empty dictionary for each of the chr
for item in range(0,22):
	chrs_len[chr_list[item]]=chr_lengths[item]

chomothriptic_chr_len = chrs_len[chromosome] # chr6 is the chromothriptic chromosome
del chrs_len[chromosome] #remove the chromothriptic chromosome 

#determine the average frequency of SVs on non-chromothriptic chromosomes - divide by 2 as there are two alleles
freq_of_svs = [] 
for chr in chrs:
	frequeny_of_SVs_tmp = chrs[chr]/chrs_len[chr]
	freq_of_svs.append(frequeny_of_SVs_tmp)

ave_sv_freq = sum(freq_of_svs)/(len(freq_of_svs)*2) #on non_chromothriptic chr (divide by 2 on)
ave_sv_freq_chromothripsis = (chomothriptic_chr-(ave_sv_freq*chomothriptic_chr_len))/chomothriptic_chr_len
#determine consistency of phasing within a block
x=[]
y=[]
z=[]
for haploblock in phasing_dict:
	#no information in the phase blocks - everything is 0
	if phasing_dict[haploblock]['Number_LOH_Block_1']==0 and phasing_dict[haploblock]['Number_LOH_Block_2']==0 and phasing_dict[haploblock]['Number_SVs_Block_1']==0 and phasing_dict[haploblock]['Number_SVs_Block_2']==0:
		phasing_dict[haploblock]['mapping']='no_info'
		y.append(haploblock)
	#no SV info but a couple of snps are present 
	if phasing_dict[haploblock]['Number_LOH_Block_1']<10 and phasing_dict[haploblock]['Number_LOH_Block_2']<10 and phasing_dict[haploblock]['Number_SVs_Block_1']==0 and phasing_dict[haploblock]['Number_SVs_Block_2']==0:
		phasing_dict[haploblock]['mapping']='no_info'
		y.append(haploblock)
	elif phasing_dict[haploblock]['Number_SVs_Block_1']/phasing_dict[haploblock]['block_length'] <= ave_sv_freq and phasing_dict[haploblock]['Number_SVs_Block_2']/phasing_dict[haploblock]['block_length'] <= ave_sv_freq:
		#if abs(phasing_dict[haploblock]['Number_SVs_Block_1'] - phasing_dict[haploblock]['Number_SVs_Block_2']) <= 1:
		phasing_dict[haploblock]['mapping']='no_info'
		y.append(haploblock)
	#chromothrpisis on haplotype 1
	elif phasing_dict[haploblock]['Number_SVs_Block_1']/phasing_dict[haploblock]['block_length'] >= ave_sv_freq:
		if phasing_dict[haploblock]['Number_SVs_Block_2']/phasing_dict[haploblock]['block_length'] <= ave_sv_freq or phasing_dict[haploblock]['Number_SVs_Block_2'] == 1:
			if phasing_dict[haploblock]['Number_SVs_Block_2']<=2:
				phasing_dict[haploblock]['mapping']='consistent_chromothripsis_1'
				x.append(haploblock)
			#if the number of SVs in one of the blocks is greater than the other with a frequency that is consistent with chromothripsis
			elif phasing_dict[haploblock]['Number_SVs_Block_1']*(ave_sv_freq/ave_sv_freq_chromothripsis) > phasing_dict[haploblock]['Number_SVs_Block_2']:
				phasing_dict[haploblock]['mapping']='consistent_chromothripsis_1'
				x.append(haploblock)
			else:
				z.append(haploblock)
		#no information because the block is too small
		elif phasing_dict[haploblock]['Number_SVs_Block_1']==1 and phasing_dict[haploblock]['Number_SVs_Block_2'] < 10:
			phasing_dict[haploblock]['mapping']='no_info'
			y.append(haploblock)
		#no information because the block is too small
		elif phasing_dict[haploblock]['Number_SVs_Block_2']==1 and phasing_dict[haploblock]['Number_SVs_Block_1'] < 10:
			phasing_dict[haploblock]['mapping']='no_info'
			y.append(haploblock)
		else:
			z.append(haploblock)
	#chromothrpisis on haplotype 2
	elif phasing_dict[haploblock]['Number_SVs_Block_2']/phasing_dict[haploblock]['block_length'] >= ave_sv_freq:
		if phasing_dict[haploblock]['Number_SVs_Block_1']/phasing_dict[haploblock]['block_length'] <= ave_sv_freq or phasing_dict[haploblock]['Number_SVs_Block_1'] == 1:
			if phasing_dict[haploblock]['Number_SVs_Block_1']<=2:
				phasing_dict[haploblock]['mapping']='consistent_chromothripsis_2'
				x.append(haploblock)
			#if the number of SVs in one of the blocks is greater than the other with a frequency that is consistent with chromothripsis
			elif phasing_dict[haploblock]['Number_SVs_Block_2']*(ave_sv_freq/ave_sv_freq_chromothripsis) > phasing_dict[haploblock]['Number_SVs_Block_1']:
				phasing_dict[haploblock]['mapping']='consistent_chromothripsis_2'
				x.append(haploblock)
			else:
				z.append(haploblock)
		#no information because the block is too small
		elif phasing_dict[haploblock]['Number_SVs_Block_1']==1 and phasing_dict[haploblock]['Number_SVs_Block_2'] < 10:
			phasing_dict[haploblock]['mapping']='no_info'
			y.append(haploblock)
		#no information because the block is too small
		elif phasing_dict[haploblock]['Number_SVs_Block_2']==1 and phasing_dict[haploblock]['Number_SVs_Block_1'] < 10:
			phasing_dict[haploblock]['mapping']='no_info'
			y.append(haploblock)
		else:
			z.append(haploblock)
	else:
		z.append(haploblock)

#for each read print the phasing for each block 
for read in bam.fetch('chr6'):
	read.qname 
	tags = read.get_tags("HP")
	chr = read.reference_name
	ref_start = read.reference_start
	ref_end = read.reference_end
	hp_tag_tmp = [item for item in tags if item[0] == "HP"]
	ps_tag_tmp = [item for item in tags if item[0] == "PS"]
	if len(hp_tag_tmp)==0:
		hp_tag = "none"
	else:
		hp_tag=hp_tag_tmp[0][1]
	if len(ps_tag_tmp)==0:
		ps_tag = "none"
	else:
		ps_tag=ps_tag_tmp[0][1]
	reads_dict[read.qname].setdefault("block",[]) #make an empty list for the blocks key as long as it doesn't already exist
	reads_dict[read.qname]["block"].append( {"chr":chr, "ref_start":ref_start, "ref_end":ref_end, "hp_tag":hp_tag, "ps_tag":ps_tag} )

#determine which read falls into category unphased, consistent_to_1, consistent_to_2 or inconsistent
for read in reads_dict:
	ps = []
	hp = []
	for segment in range(len(reads_dict[read]['block'])):
		hp.append(reads_dict[read]['block'][segment]['hp_tag'])
		ps.append(reads_dict[read]['block'][segment]['ps_tag'])
	if all(phase == "none" for phase in hp):
		reads_dict[read]['haplotype']="unphased"
	elif all(phase == 1 or phase == "none" for phase in hp):
		reads_dict[read]['haplotype']="consistent_to_1"
	elif all(phase == 2 or phase == "none" for phase in hp):
		reads_dict[read]['haplotype']="consistent_to_2"
	elif 1 in hp and 2 in hp:
		reads_dict[read]['haplotype']="inconsistent_phasing"
		print('dont_know')
	else:
		read
	reads_dict[read]['phase_set'] = list(set(ps))		
	reads_dict[read]['phase_set'] = list(filter(lambda a: a != "none", reads_dict[read]['phase_set']))

#QC check how many reads fall into each catagory
unphased = 0
consistent_to_1 = 0
consistent_to_2 = 0
inconsistent_phasing = 0
for read in reads_dict:
	if reads_dict[read]['haplotype']=="unphased":
		unphased=unphased+1
	if reads_dict[read]['haplotype']=="consistent_to_1":
		consistent_to_1=consistent_to_1+1
	if reads_dict[read]['haplotype']=="consistent_to_2":
		consistent_to_2 =consistent_to_2+1
	if reads_dict[read]['haplotype']=="inconsistent_phasing":
		inconsistent_phasing=inconsistent_phasing+1

#find the reads which span multiple phase sets 
reads_mapping_multiple_phase_sets=[]
for read in reads_dict:
	if len(reads_dict[read]['phase_set'])>1:
		reads_mapping_multiple_phase_sets.append(read)


missing_blocks = [] # need to look at these
chromothriptic_reads = []
wild_type_reads = []
unphased_reads = [] # we will randomly assign these
inconsistent_reads = [] # should be none 
no_info = [] # split by block and then randomly assign 
for read in reads_dict:
	if reads_dict[read]['haplotype'] == "unphased":
		unphased_reads.append(read)
	elif reads_dict[read]['haplotype'] == "inconsistent_phasing": #none of them
		inconsistent_reads.append(read)
	elif reads_dict[read]['haplotype'] == "consistent_to_1" or  reads_dict[read]['haplotype'] == "consistent_to_2":
		haploblock = reads_dict[read]['phase_set'][0]
		if reads_dict[read]['phase_set'][0] in phasing_dict:
			if "mapping" in phasing_dict[haploblock].keys(): 
				if phasing_dict[haploblock]['mapping']=="no_info":
					#do something to split reads
					no_info.append(read)
				elif reads_dict[read]['haplotype'] == "consistent_to_1":
					if phasing_dict[haploblock]['mapping'] =='consistent_chromothripsis_1':
						chromothriptic_reads.append(read)
					if phasing_dict[haploblock]['mapping'] =='consistent_chromothripsis_2':
						wild_type_reads.append(read)
				elif reads_dict[read]['haplotype'] == "consistent_to_2":
					if phasing_dict[haploblock]['mapping'] =='consistent_chromothripsis_2':
						chromothriptic_reads.append(read)
					if phasing_dict[haploblock]['mapping'] =='consistent_chromothripsis_1':
						wild_type_reads.append(read)
			else:
				missing_blocks.append(haploblock)
		else:
			phasing_dict[haploblock]={'Number_SNPs':0, 'Number_LOH_Block_1':0, 'Number_LOH_Block_2':0, 'Number_SVs_Block_1':0, 'Number_SVs_Block_2':0, 'mapping':'no_info'}
			no_info.append(read)
	else:
		read

#look at the missing blocks since they contain phasing errors
missing_blocks = list(set(missing_blocks))

#list of LOH regions  
loh_regions_input = pd.read_csv(loh_regions_file, sep="\t", skiprows=0, header=None) 
loh_regions_input.columns = ['chr','start','end']
#make a dictionary so the data is more accessible
loh_regions_dict={x:{} for x in list(loh_regions_input['chr'])} #make an empty dictionary for each chromosome
for key in loh_regions_dict.keys():
	loh_regions_dict[key]={"start":[], "end":[]}

for key in range(len(loh_regions_input)):
	loh_regions_dict[loh_regions_input['chr'][key]]['start'].append(loh_regions_input['start'][key])
	loh_regions_dict[loh_regions_input['chr'][key]]['end'].append(loh_regions_input['end'][key])


#get the names of the blocks which have no info
no_info_blocks=[]
for read in no_info:
	if reads_dict[read]['phase_set'][0] not in no_info_blocks:
		no_info_blocks.append(reads_dict[read]['phase_set'][0])

#set up a dictions which you can fill with read names
reads_by_block = { x:{} for x in no_info_blocks }
for block in reads_by_block:
	reads_by_block[block]["haplotype1"]=[]
	reads_by_block[block]["haplotype2"]=[]

#fill the dictionary with read names
for read in no_info:
	ps = reads_dict[read]['phase_set'][0]
	if reads_dict[read]['haplotype'] == "consistent_to_1":
		reads_by_block[ps]["haplotype1"].append(read)
	if reads_dict[read]['haplotype'] == "consistent_to_2":
		reads_by_block[ps]["haplotype2"].append(read)

x = 0
y = 0 
z = 0 
a = 0
b = 0
chromothriptic_no_info=[]
wild_type_no_info=[]
for block in reads_by_block:
	haplotype1_sv_reads = 0
	haplotype2_sv_reads = 0
	for read in reads_by_block[block]["haplotype1"]:
		if read in reads_with_SVs:
			haplotype1_sv_reads=haplotype1_sv_reads+1
	for read in reads_by_block[block]["haplotype2"]:
		if read in reads_with_SVs:
			haplotype1_sv_reads=haplotype2_sv_reads+1
	haplotype1_sv_reads
	haplotype2_sv_reads
	print("---")
	#deal with SVs
	if haplotype2_sv_reads*3 < haplotype1_sv_reads  and haplotype2_sv_reads!=1 and haplotype1_sv_reads!=1:
		chromothriptic_no_info.append(reads_by_block[block]['haplotype1'])
		wild_type_no_info.append(reads_by_block[block]['haplotype2'])
		z=z+1
		print("z")
	#deal with SVs
	elif haplotype1_sv_reads*3 < haplotype2_sv_reads and haplotype1_sv_reads!=1 and haplotype1_sv_reads!=1:
		chromothriptic_no_info.append(reads_by_block[block]['haplotype2'])
		wild_type_no_info.append(reads_by_block[block]['haplotype1'])
		a=a+1
		print("a")
	#deal with LOH
	elif phasing_dict[block]['Number_LOH_Block_1']==0 and phasing_dict[block]['Number_LOH_Block_2']>1:
		chromothriptic_no_info.append(reads_by_block[block]['haplotype1'])
		wild_type_no_info.append(reads_by_block[block]['haplotype2'])
		x=x+1
		print("x")
	elif phasing_dict[block]['Number_LOH_Block_2']==0 and phasing_dict[block]['Number_LOH_Block_1']>1:
		chromothriptic_no_info.append(reads_by_block[block]['haplotype2'])
		wild_type_no_info.append(reads_by_block[block]['haplotype1'])
		y=y+1
		print("y")
	#randomly assign everything else
	else:
		b=b+1
		print("b")
		chromothriptic_block = random.choice(list(reads_by_block[block].keys()))
		if chromothriptic_block == 'haplotype1':
			chromothriptic_no_info.append(reads_by_block[block]['haplotype1'])
			wild_type_no_info.append(reads_by_block[block]['haplotype2'])
		else:
			chromothriptic_no_info.append(reads_by_block[block]['haplotype2'])
			wild_type_no_info.append(reads_by_block[block]['haplotype1'])

chromothriptic_no_info_list = [y for x in chromothriptic_no_info for y in x]
wild_type_no_info_list = [y for x in wild_type_no_info for y in x]

#determin if the unphased reads have SV
chromothriptic_unphased=[]
wild_type_unphased=[]
unphased_reads_without_SVs=[]
for read in unphased_reads:
	if read in reads_with_SVs:
		chromothriptic_unphased.append(read)
	else:
		unphased_reads_without_SVs.append(read)

#determine how many reads shoud be in each file and then shuffle the readnames and split then randomly 
split_for_chromothriptic=round((len(unphased_reads)/2)-len(chromothriptic_unphased))
shuffled_unphased_reads = unphased_reads.copy()
random.shuffle(shuffled_unphased_reads)
for read in range(len(shuffled_unphased_reads[:split_for_chromothriptic])):
	chromothriptic_unphased.append(unphased_reads_without_SVs[read])

for read in range(split_for_chromothriptic, len(unphased_reads_without_SVs)):
	wild_type_unphased.append(unphased_reads_without_SVs[read])

chromothriptic_reads_final = chromothriptic_reads + chromothriptic_no_info_list + chromothriptic_unphased
wildtype_reads_final = wild_type_reads + wild_type_no_info_list + wild_type_unphased

missing_blocks = [] # need to look at these

#output reads to two files:
obam = pysam.AlignmentFile("{}{}{}{}".format(outputdir,"/chromothriptic_reads_all_",chromosome,"_ccs.sam"), "w", template=bam)

for b in bam.fetch(chromosome):
    if b.query_name in chromothriptic_reads_final:
        obam.write(b)

obam.close()

outbam = pysam.AlignmentFile("{}{}{}{}".format(outputdir,"/wildtype_reads_all_",chromosome,"_ccs.sam"), "w", template=bam)
for b in bam.fetch(chromosome):
    if b.query_name in wildtype_reads_final:
        outbam.write(b)

outbam.close()
bam.close()
