# /home/users/pjh/scripts/annotation/SV/gremlin_step2/02_SV_info_edit.v3.py

# args[1]: sv vcf  # must be sorted. (del 3-5, dup 5-3, inv pos1<pos2, tra chr1<chr2; 1...22,X,Y,MT)
# args[2]: tumor bam
# args[3]: normal bam
# args[4]: reference.fasta.fai

#import sys
#import os
#import collections
#import itertools
#import re
#import random
#import math
import argparse
#from statistics import median

#import pyranges as pr
import pysam

import sv_functions as sv

def main():
	
	# argument parsing
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', help = 'SV vcf file')
	parser.add_argument('tbam', help = 'tumor bam file')
	parser.add_argument('nbam', help = 'normal bam file')
	parser.add_argument('fai', help = 'reference fai file')
	parser.add_argument('refver', help = 'reference genome version, must be 19 or 38')
	parser.add_argument('--include_nonprimary', action = 'store_true', help = 'If set, non-primary contigs are included in analysis')
	args = parser.parse_args()

	in_file = open(args.infile)
	tbam_file = pysam.AlignmentFile(args.tbam,'rb')
	nbam_file = pysam.AlignmentFile(args.nbam,'rb')
	out_file = open(args.infile + '.SVinfo', mode = 'w', buffering = 1)
	fai_path = args.fai
	refver = args.refver
	nonprim = args.include_nonprimary # default is Fault, not to include non-primary contigs

	# get chromosome size dictionary
	chr_size, chr_list = sv.get_chr_size(fai_path)

	# set constant parameters
	satellite_path,\
	pr_satellite,\
	primary_contigs,\
	r_limit,\
	r_limit_dp,\
	terinfo_list = sv.set_constants(refver)
	
	# write header line to the output file
	firstline = next(in_file)
	assert firstline.startswith('CHR1')
	hdr = sv.get_header(firstline)
	out_file.write('\t'.join(hdr) + '\n')
	
	# step 2-7
	for line in in_file:
		sv.printErr(line, end='')

		linesp = line.strip().split('\t')
		chr1, pos1, chr2, pos2, sv_type, terinfo = linesp
		pos1 = int(pos1) ; pos2 = int(pos2)
		dist = pos2 - pos1
		ter1, ter2 = terinfo.split('to')
		
		# variant filtering
		if sv.variant_filter(chr1, pos1, chr2, pos2, pr_satellite, primary_contigs, sv_type, nonprim, chr_size): # Filter out if True is returned
			continue
		
		# sanity check
		if sv_type!='TRA' and dist < 0:
			sv.printErr('Sorting error')
			sv.printErr(line, end='') ; sys.exit(1)
		if terinfo not in terinfo_list:
			sv.printErr('6th column entry is not one of "3to3", "3to5", "5to3", or "5to5".')
			sv.printErr(line, end='') ; sys.exit(1)
		if sv_type == 'INV' and not ((ter1=='3' and ter2=='3') or (ter1=='5' and ter2=='5')):
			sv.printErr('SV type is INV but 6th column is not "3to3" nor "5to5".')
			sv.printErr(line, end='') ; sys.exit(1)
			
		# get rplists
		rplist_dict = \
			sv.get_rplist(
				sv.collect_fetchrange(chr1, pos1, ter1, chr2, pos2, ter2, sv_type, dist, chr_size),
				tbam_file, nbam_file, ter1, r_limit_dp, chr_list, primary_contigs, nonprim
			)
			
		#(1) BP information
		t_info, n_info = sv.step1(chr1, pos1, ter1, chr2, pos2, ter2, sv_type, dist, chr_list, rplist_dict, r_limit)
		
		#(2) BP adjustment
		chr1, pos1, chr2, pos2, info_list = sv.step2(chr1, pos1, chr2, pos2, terinfo, sv_type, t_info, n_info, line)
		
		# run calcfinalcount (find_discordant_reads) - result of this object will be used also in step 6
		wrapper_calcFinalCount = sv.calc_final_count(chr1=chr1, pos1=pos1, ter1=ter1, chr2=chr2, pos2=pos2, ter2=ter2,
													rplist_dict=rplist_dict, chr_list=chr_list,
													r_limit=r_limit)
		wrapper_calcFinalCount.main()
			
		#(3) SV vaf
		print_list = sv.step3(chr1, pos1, ter1, chr2, pos2, ter2, rplist_dict, wrapper_calcFinalCount, r_limit)
		
		#(4) Paired normal soft clip
		pnsc = sv.step4(chr1, pos1, ter1, chr2, pos2, ter2, rplist_dict, r_limit)
		
		# (5) Background information
		tres1, tres2, nres1, nres2 = sv.step5(chr1, pos1, chr2, pos2, rplist_dict, r_limit)
		
		# (6) Mapping quality and position variation
		mqpos_list = sv.step6(wrapper_calcFinalCount)
		
		# (7) Depth ratio change
		depth_ratio_change_bp1, depth_ratio_change_bp2 = sv.step7(chr1, pos1, ter1, chr2, pos2, ter2, tbam_file, nbam_file, chr_size, searchlen = 500)
		
		# write results to output file
		out_file.write(
			'\t'.join([
				line.strip(),
				','.join(t_info) if t_info != None else 'NA',
				','.join(n_info) if n_info != None else 'NA',
				'\t'.join(info_list),
				'\t'.join(print_list),
				str(pnsc),
				tres1, 
				tres2, 
				nres1, 
				nres2,
				'\t'.join(mqpos_list),
				str(depth_ratio_change_bp1),
				str(depth_ratio_change_bp2)
			]) + '\n'
		)
		
if __name__ == '__main__':
	main()
