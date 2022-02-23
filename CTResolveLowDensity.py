#! /usr/bin/python

while read CHR; do
    bsub -G team154-grp -oo logs/${CHR}_phase_CCS_reads.out -q long -R'select[mem>12000] rusage[mem=12000]' -M12000 -n1 -J ${CHR}_phase_CCS_reads \
    python CTResolveLowDensity.py --chromosome ${CHR} \
    --bam /lustre/scratch117/casm/team154/ji2/OesophagealLongRead_WTSI-OESO_117/phaseCCS/whatshap/phased_bams/haplotag_${CHR}_phased_ccs_no_simple_repeats.bam \
    --readnames ../OesophagealLongRead_WTSI-OESO_117/phaseCCS/phased_CCS_reads/allreadnames/allreadnames_${CHR}.txt \
    --haplotype1_reads ../OesophagealLongRead_WTSI-OESO_117/phaseCCS/whatshap/phased_reads/OESO_117_hg38_${CHR}_tagged_reads_haplotype1_readnames_chr_pos.txt \
    --haplotype2_reads ../OesophagealLongRead_WTSI-OESO_117/phaseCCS/whatshap/phased_reads/OESO_117_hg38_${CHR}_tagged_reads_haplotype2_readnames_chr_pos.txt \
    --haplotype_blocks  ../OesophagealLongRead_WTSI-OESO_117/phaseCCS/whatshap/whatshap_stats/OESO_117_${CHR}_ccs_blocks.tsv \
    --phased_haplotype_blocks /lustre/scratch117/casm/team154/pc8/Jannat_pacbio/OESO_117/OESO_117_haplotype_blocks.txt  \
    --loh_regions tmp \
	--outputdir phased_CCS_reads/
done < chromosomes.txt




import pysam 
import pandas as pd
import random
import os 
from os import system
import argparse 

parser = argparse.ArgumentParser(description='usage:  CTResolveLowDensity.py --chromosome chromosome --bam bam --readnames readnames\
	--haplotype1_reads haplotype1_reads --haplotype2_reads haplotype2_reads --haplotype_blocks haplotype_blocks --phased_haplotype_blocks phased_haplotype_blocks \
	--loh_regions --loh_regions --outputdir path')
parser.add_argument('--chromosome', help='chromosome',required=True)
parser.add_argument('--bam', help='path to bamfile',required=True)
parser.add_argument('--readnames', help='path to readnames list',required=True)
parser.add_argument('--haplotype1_reads', help='readnames pos chr haplotype1 ccs',required=True)
parser.add_argument('--haplotype2_reads', help='readnamespos chr  haplotype2 ccs',required=True)
parser.add_argument('--haplotype_blocks', help='path to the haplotype blocks produced by whatshap stats',required=True)
parser.add_argument('--phased_haplotype_blocks', help='path to the phased haplotype blocks produced by Peters code',required=True)
parser.add_argument('--loh_regions', help='path to loh file',required=True)
parser.add_argument('--outputdir', help='path to the output directory',required=True)
args = parser.parse_args()


chromosome=args.chromosome
bamfile=args.bam
readnames_file=args.readnames
haplotype1_readnames_chr_pos_file=args.haplotype1_reads
haplotype2_readnames_chr_pos_file=args.haplotype2_reads
haplotype_blocks_file=args.haplotype_blocks
phasing_haplotype_blocksfile=args.phased_haplotype_blocks
loh_regions_file=args.loh_regions
outputdir=args.outputdir


##LOAD DATA
#bam
bam = pysam.AlignmentFile(bamfile, "rb") #load bam file


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
phasing = pd.read_csv(phasing_haplotype_blocksfile, sep="\t", skiprows=1, header=None)
phasing.columns = ["chr","Haplotype","Number_SNPs","ave_VAF_block1","ave_depth_block1","ave_depth_block2"]
#make a dictionary so the data is more accessible
phasing_dict={} #make an empty dictionary for each of the SVs
for key in range(len(phasing)):
    if phasing['chr'][key]==chromosome:
        phasing_dict[phasing['Haplotype'][key]]={}

for key in range(len(phasing)):
	if phasing['chr'][key]==chromosome:
		phasing_dict[phasing['Haplotype'][key]]={'Number_SNPs':phasing['Number_SNPs'][key], 'ave_VAF_block1':phasing['ave_VAF_block1'][key], 'ave_VAF_block2':1-phasing['ave_VAF_block1'][key],  'ave_depth_block1':phasing['ave_depth_block1'][key], 'ave_depth_block2':phasing['ave_depth_block2'][key]}


#find the block length, start and end and then match them with the phasing info
for key in phasing_dict:
	for haplotype in haplotype_blocks_dict:
		if key == str(haplotype):
			phasing_dict[key]['block_length']=haplotype_blocks_dict[int(key)]['length']
			phasing_dict[key]['from']=haplotype_blocks_dict[int(key)]['from']
			phasing_dict[key]['to']=haplotype_blocks_dict[int(key)]['to']

#ave_sv_freq = sum(freq_of_svs)/(len(freq_of_svs)*2) #on non_chromothriptic chr (divide by 2 on)
#ave_sv_freq_chromothripsis = (chomothriptic_chr-(ave_sv_freq*chomothriptic_chr_len))/chomothriptic_chr_len
#determine consistency of phasing within a block
x=[]
y=[]
z=[]
for haploblock in phasing_dict:
	#no information in the phase blocks - everything is 0
	if phasing_dict[haploblock]['ave_VAF_block1']==0 and phasing_dict[haploblock]['ave_VAF_block2']==1 and phasing_dict[haploblock]['ave_depth_block1']==0 and phasing_dict[haploblock]['ave_depth_block2']==0:
		phasing_dict[haploblock]['mapping']='no_info'
		y.append(haploblock)
	#no diffence in VAF info so we randomly assign keep the block in haplotype1 
	elif 0.45<=phasing_dict[haploblock]['ave_VAF_block1']<=0.55:
		phasing_dict[haploblock]['mapping']='no_info'
		y.append(haploblock)
	#rearrangments on haplotype 1 if the copy number is over 2 then rearranged, if it is LOH also rearranged
	elif 0.55<phasing_dict[haploblock]['ave_VAF_block1']<0.9:
		phasing_dict[haploblock]['mapping']='consistent_rearrangement_1'
		x.append(haploblock)
	elif phasing_dict[haploblock]['ave_VAF_block1']<0.1:
		phasing_dict[haploblock]['mapping']='consistent_rearrangement_1'
		x.append(haploblock)
	#rearrangments on haplotype 2 if the copy number is over 2 then rearranged, if it is LOH also rearranged
	elif 0.55<phasing_dict[haploblock]['ave_VAF_block2']<0.9:
		phasing_dict[haploblock]['mapping']='consistent_rearrangement_2'
		x.append(haploblock)
	elif phasing_dict[haploblock]['ave_VAF_block2']<0.1:
		phasing_dict[haploblock]['mapping']='consistent_rearrangement_2'
		x.append(haploblock)
	else:
		z.append(haploblock)


print(len(x))
print(len(y))
print(len(z))

#for each read print the phasing for each block 
for read in bam.fetch(chromosome):
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
haplotype1_reads = []
haplotype2_reads = []
unphased_reads = [] # we will randomly assign these
inconsistent_reads = [] # should be none 
no_info = [] # split by block and then randomly assign 
for read in reads_dict:
	if reads_dict[read]['haplotype'] == "unphased":
		unphased_reads.append(read)
	elif reads_dict[read]['haplotype'] == "inconsistent_phasing": #none of them
		inconsistent_reads.append(read)
	elif reads_dict[read]['haplotype'] == "consistent_to_1" or  reads_dict[read]['haplotype'] == "consistent_to_2":
		haploblock = str(reads_dict[read]['phase_set'][0])
		if haploblock in phasing_dict:
			if "mapping" in phasing_dict[haploblock].keys(): 
				if phasing_dict[haploblock]['mapping']=="no_info":
					#do something to split reads
					no_info.append(read)
				elif reads_dict[read]['haplotype'] == "consistent_to_1":
					if phasing_dict[haploblock]['mapping'] =='consistent_rearrangement_1':
						haplotype1_reads.append(read)
					if phasing_dict[haploblock]['mapping'] =='consistent_rearrangement_2':
						haplotype2_reads.append(read)
				elif reads_dict[read]['haplotype'] == "consistent_to_2":
					if phasing_dict[haploblock]['mapping'] =='consistent_rearrangement_2':
						haplotype1_reads.append(read)
					if phasing_dict[haploblock]['mapping'] =='consistent_rearrangement_1':
						haplotype2_reads.append(read)
			else:
				missing_blocks.append(haploblock)
		else:
			phasing_dict[haploblock]={'Number_SNPs':0, 'ave_VAF_block1':0, 'ave_VAF_block2':0, 'ave_depth_block1':0, 'ave_depth_block2':0, 'mapping':'no_info'}
			no_info.append(read)
	else:
		read

#look at the missing blocks since they contain phasing errors
missing_blocks = list(set(missing_blocks))

print("number of inconsistent blocks = ",inconsistent_phasing )

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

#set up a dictionary which you can fill with read names
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


b = 0
haplotype1_no_info=[]
haplotype2_no_info=[]
for block in reads_by_block:
	b=b+1
	its_loh='no'
	if chromosome in loh_regions_dict:
		for segment in range(len(loh_regions_dict[chromosome]['start'])):
			if loh_regions_dict[chromosome]['start'][segment]<block<loh_regions_dict[chromosome]['end'][segment]:
				its_loh='yes'
	if its_loh=='yes':
		haplotype2_reads=haplotype2_reads+reads_by_block[block]['haplotype1']+reads_by_block[block]['haplotype2']
	else:
		haplotype1_block = random.choice(list(reads_by_block[block].keys()))
		if haplotype1_block == 'haplotype1':
			haplotype1_no_info.append(reads_by_block[block]['haplotype1'])
			haplotype2_no_info.append(reads_by_block[block]['haplotype2'])
		else:
			haplotype1_no_info.append(reads_by_block[block]['haplotype2'])
			haplotype2_no_info.append(reads_by_block[block]['haplotype1'])

haplotype1_no_info_list = [y for x in haplotype1_no_info for y in x]
haplotype2_no_info_list = [y for x in haplotype2_no_info for y in x]

haplotype2_loh=[]
haplotype1_unphased=[]
haplotype2_unphased=[]
#determine how many reads shoud be in each file and then shuffle the readnames and split then randomly 
reads_to_remove=[]
for read in unphased_reads:
	if chromosome in loh_regions_dict:
		for segment in range(len(loh_regions_dict[chromosome]['start'])):
			for block in range(len(reads_dict[read]['block'])):
				if reads_dict[read]['block'][block]['chr']==chromosome and \
				loh_regions_dict[chromosome]['start'][segment]<reads_dict[read]['block'][block]['ref_start']<loh_regions_dict[chromosome]['end'][segment] and \
				loh_regions_dict[chromosome]['start'][segment]<reads_dict[read]['block'][block]['ref_end']<loh_regions_dict[chromosome]['end'][segment]:
					haplotype2_loh.append(read)
					reads_to_remove.append(read)

reads_to_remove=set(reads_to_remove)
for read in reads_to_remove:
	unphased_reads.remove(read)

split=round(len(unphased_reads)/2)
shuffled_unphased_reads = unphased_reads.copy()
random.shuffle(shuffled_unphased_reads)
for read in range(len(shuffled_unphased_reads[:split])):
	haplotype1_unphased.append(shuffled_unphased_reads[read])

for read in range(split, len(unphased_reads)):
	haplotype2_unphased.append(shuffled_unphased_reads[read])


haplotype1_reads_final = haplotype1_reads + haplotype1_no_info_list + haplotype1_unphased
haplotype2_reads_final = haplotype2_reads + haplotype2_no_info_list + haplotype2_unphased + haplotype2_loh

missing_blocks = [] # need to look at these

#output reads to two files:
obam = pysam.AlignmentFile("{}{}{}{}".format(outputdir,"/haplotype1_reads_all_",chromosome,"_ccs.sam"), "w", template=bam)

for b in bam.fetch(chromosome):
    if b.query_name in haplotype1_reads_final:
        obam.write(b)

obam.close()

outbam = pysam.AlignmentFile("{}{}{}{}".format(outputdir,"/haplotype2_reads_all_",chromosome,"_ccs.sam"), "w", template=bam)
for b in bam.fetch(chromosome):
    if b.query_name in haplotype2_reads_final:
        outbam.write(b)

outbam.close()
bam.close()

