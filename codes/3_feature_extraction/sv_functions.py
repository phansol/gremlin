# /home/users/pjh/scripts/annotation/SV/gremlin_step2/sv_functions_edit.py
# /home/users/pjh/scripts/annotation/SV/gremlin_step2/02_SV_info_edit.v3.py

import sys
import os
import collections
import itertools
import re
import random
import math
from statistics import median
import pyranges as pr
import pysam

cigarpat = re.compile('([0-9]+)([MIDNSHP=XB])')
cigarOpDict = dict(zip(('M','I','D','N','S','H','P','=','X','B'), range(10)))
SA_dict_keys = ['chr', 'pos', 'strand', 'cigarstring', 'MQ', 'NM']
readlen_cigarops = [0, 1, 4, 5]

def make_cigartuple(cigarstring):
	return [ (cigarOpDict[x[1]], int(x[0])) for x in cigarpat.findall(cigarstring) ]

def change_chr_to_int(chr1):
	chr1=chr1.replace('chr','')
	chr_n=int(chr1.replace('X','23').replace('Y','24'))
	return(chr_n)

# defining readplus class
class readplus:
	def __init__(self, read, ter1, chr_list, primary_contigs, nonprim):
		self.sc_co = 5
		
		self.read = read
		self.ter1 = ter1

		#self.chr_n = change_chr_to_int(self.read.reference_name)
		self.chr_n = chr_list.index(self.read.reference_name)

		self.read_size = self.read.infer_read_length() # this includes hard-clip lengths
		self.M_range = self.find_M_range(self.read.cigartuples)
		self.M_range_len = self.M_range[1] - self.M_range[0]
		self.get_SA_list(chr_list, primary_contigs, nonprim)
		
		if self.ter1 == '5':
			self.sc_seq = self.read.query_sequence[self.read.cigartuples[0][1]-1-self.sc_co+1 : self.read.cigartuples[0][1]-1+1]
		elif self.ter1 == '3':
			self.sc_seq = self.read.query_sequence[self.read.cigartuples[-1][1]*(-1) : self.read.cigartuples[-1][1]*(-1)+self.sc_co]
			
		
	def get_SA_list(self, chr_list, primary_contigs, nonprim): # entries of SA_list: dict with SA_dict_keys as keys
		self.SA_list = list()
		if self.read.has_tag('SA'):
			for SA_indi in self.read.get_tag('SA').split(';')[:-1]:
				SAdic = dict(zip(SA_dict_keys, SA_indi.split(',')))

				# skip if SA chromosome is not a primary contig
				if not nonprim:
					if SAdic['chr'] not in primary_contigs: # added on 210128
						continue

				SAdic['pos'] = int(SAdic['pos'])
				SAdic['MQ'] = int(SAdic['MQ'])
				SAdic['NM'] = int(SAdic['NM'])
				
				#SAdic['chr_n'] = change_chr_to_int(SAdic['chr'])
				SAdic['chr_n'] = chr_list.index(SAdic['chr'])

				SAdic['is_reverse'] = True if SAdic['strand'] == '-' else False
				SAdic['cigartuples'] = make_cigartuple(SAdic['cigarstring'])
				SAdic['M_range'] = self.find_M_range(SAdic['cigartuples'])
				SAdic['M_range_len'] = SAdic['M_range'][1] - SAdic['M_range'][0]
				SAdic['current_m'] = sum([ x[1] for x in SAdic['cigartuples'] if x[0] == 0 ])
				SAdic['current_d'] = sum([ x[1] for x in SAdic['cigartuples'] if x[0] == 2 ])
				SAdic['read_size'] = sum([ x[1] for x in SAdic['cigartuples'] if x[0] in readlen_cigarops ])
				
				self.SA_list.append(SAdic)
			
	def find_M_range(self, cigartuples):
		m_start = m_end = 0  # m_start: just before the start, m_end= the exact end
		m_count=0
		for (t, n) in cigartuples:
			if t == 0:
				m_count +=1
		if m_count ==1:
			for (t,n) in cigartuples:
				if t!=0 and t!=1:
					m_start+=n
				elif t==0:
					m_end=m_start+n
					break
		elif m_count > 1:
			find_m=0;m_length=0
			for (t,n) in cigartuples:
				if find_m==0 and t!=0 and t!=1:
					m_start+=n
				elif find_m >0 and t!=0 and t!=1:
					m_length+=n
				elif t==0:
					find_m+=1
					if find_m < m_count:
						m_length+=n
					elif find_m == m_count:
						m_end=m_start+m_length
						break
		return (m_start, m_end)
			
			
# functions used in multiple classes
def find_discordant_reads_poscalc(chr1, pos1, ter1, chr2, pos2, ter2, fors, bacs):
	if ter1=='3':
		pos1_start = pos1 - fors ; pos1_end = pos1 + bacs
	elif ter1=='5': 
		pos1_start = pos1 - bacs ; pos1_end = pos1 + fors
	if ter2=='3': 
		pos2_start = pos2 - fors ; pos2_end = pos2 + bacs
	elif ter2=='5': 
		pos2_start = pos2 - bacs ; pos2_end = pos2 + fors
	pos1_start = max(pos1_start, 1); pos2_start = max(pos2_start, 1)
	
	if chr1 == chr2 and ter1 == '5' and ter2 == '3' and pos1 < pos2:  # exceptional short duplication
		pos1_end = min(pos1_end, pos1+(pos2-pos1)/2)
		pos2_start = max(pos2_start, pos2-(pos2-pos1)/2)
	elif chr1 == chr2 and ter1 == '3' and ter2 =='5' and pos2 < pos1:
		pos2_end = min(pos2_end, pos2+(pos1-pos2)/2)
		pos1_start = max(pos1_start, pos1-(pos1-pos2)/2)
		
	return pos1_start, pos1_end, pos2_start, pos2_end


def readfilter_1(SAdic, chr2, pos2):
	if \
	( \
		( SAdic['cigartuples'][0][0]==4 or SAdic['cigartuples'][0][0]==5 ) and \
		SAdic['chr'] == chr2 and \
		abs(SAdic['pos']-pos2) <=1 \
	) or \
	( \
		( SAdic['cigartuples'][-1][0]==4 or SAdic['cigartuples'][-1][0]==5 ) and \
		SAdic['chr'] == chr2 and \
		abs(SAdic['pos'] + SAdic['current_m'] + SAdic['current_d']-1-pos2) <=1 \
	):
		return True
	else:
		return False
	
	
def readpair_to_SV(
	chr1, pos1, isrev1, M_range1, M_range_len1, read_size1, chr_n1, 
	chr2, pos2, isrev2, M_range2, M_range_len2, read_size2, chr_n2
):
	# determine whether the two have same directions
	samedir = isrev1 == isrev2
	
	# modify M_range1 if different directions
	if not samedir:
		M_range1 = [read_size1 - M_range1[1], read_size1 - M_range1[0]]
	
	# set default overlap value
	overlap = False
	
	if samedir:
		if \
		(M_range1[0] <= M_range2[0] and M_range1[1] >= M_range2[1]) or \
		(M_range2[0] <= M_range1[0] and M_range2[1] >= M_range1[1]):
			overlap = True
		
		else:
			if M_range1[1] > M_range2[1]:
				MHLEN=M_range2[1]-M_range1[0]
				bp1 = pos1
				bp2 = pos2 + M_range_len2 - 1
				terminal1, terminal2 = '5', '3'
				if chr1!=chr2:
					rearr="TRA"
					if chr_n1 < chr_n2: 
						ori='rs'
					elif chr_n1 > chr_n2: 
						ori='sr'
				else:
					if bp1<bp2: 
						rearr = "DUP"
						ori = 'rs'
					elif bp1>bp2: 
						rearr = "DEL"
						ori = 'sr'
					else:
						overlap = True

			elif M_range2[1] > M_range1[1]:
				MHLEN = M_range1[1]-M_range2[0]
				bp1=pos1 + M_range_len1 - 1
				bp2=pos2
				terminal1, terminal2 = '3', '5'
				if chr1 != chr2:
					rearr="TRA"
					if chr_n1 < chr_n2: 
						ori='rs'
					elif chr_n1 > chr_n2: 
						ori='sr'
				else:
					if bp1<bp2: 
						rearr="DEL"
						ori='rs'
					elif bp1>bp2: 
						rearr="DUP"
						ori='sr'
					else:
						overlap = True
						
	else:  # opposite direction
		if \
		(M_range1[0] <= M_range2[0] and M_range1[1] >= M_range2[1]) or \
		(M_range2[0] <= M_range1[0] and M_range2[1] >= M_range1[1]): 
			overlap = True
		
		else:
			if M_range1[1] > M_range2[1]:
				MHLEN = M_range2[1]-M_range1[0]
				bp1 = pos1 + M_range_len1 - 1
				bp2 = pos2 + M_range_len2 - 1
				terminal1, terminal2 = '3', '3'
				if chr1 != chr2:
					rearr="TRA"
					if chr_n1 < chr_n2: 
						ori='rs'
					elif chr_n1 > chr_n2: 
						ori='sr'
				else:
					rearr="INV"
					if bp1 < bp2: 
						ori='rs'
					elif bp1 > bp2:
						ori='sr'
					else:
						overlap = True
						
			elif M_range2[1] > M_range1[1]:
				MHLEN = M_range1[1]-M_range2[0]
				bp1 = pos1
				bp2 = pos2
				terminal1, terminal2 = '5', '5'
				if chr1 != chr2:
					rearr="TRA"
					if chr_n1 < chr_n2: 
						ori='rs'
					elif chr_n1 > chr_n2: 
						ori='sr'
				else:
					rearr="INV"
					if bp1 < bp2: 
						ori="rs"
					elif bp1 > bp2:
						ori="sr"
					else:
						overlap = True
						
	if overlap:
		return None
	else:
		return MHLEN, bp1, bp2, terminal1, terminal2, rearr, ori
	
	
	
# base class for each main job step class
class worker:
	def __init__(self, **kwargs):
		for k,v in kwargs.items():
			setattr(self, k, v)
		if 'r_limit' in dir(self):
			self.readlimit = False if self.r_limit == None else True
		
	def check_within_fetchrange(self, read, fetchrange):
		if (
			read.reference_name == fetchrange[0] and
			read.reference_start < fetchrange[2] and
			read.reference_end > fetchrange[1]
		):
			return True
		else:
			return False
		
# STEP 01 #
class step1_wrapper(worker): 
	# args: chr1, pos1, ter1, chr2, pos2, ter2, sv_type, dist, chr_list, rplist_dict, r_limit
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.fors = 700 ; self.bacs = 100
		
		
	def main(self):
		self.get_fetchrange()
		start1, end1, start2, end2, target_start1, target_end1, target_start2, target_end2 = self.fetchrange_dict['allresults']
		
		t_list1 = self.find_SA_reads(self.chr1, start1, end1, self.chr2, target_start2, target_end2, self.rplist_dict['tumor']['chr1'], self.chr_list)
		t_list2 = self.find_SA_reads(self.chr2, start2, end2, self.chr1, target_start1, target_end1, self.rplist_dict['tumor']['chr2'], self.chr_list)
		t_list = t_list1 + t_list2
		if len(t_list) > 0:
			t_dic = collections.Counter(t_list)
			t_info = [ f'{k}({v})' for (k,v) in t_dic.items() ]
		else:
			t_info = None
			
		n_list1 = self.find_SA_reads(self.chr1, start1, end1, self.chr2, target_start2, target_end2, self.rplist_dict['normal']['chr1'], self.chr_list)
		n_list2 = self.find_SA_reads(self.chr2, start2, end2, self.chr1, target_start1, target_end1, self.rplist_dict['normal']['chr2'], self.chr_list)
		n_list = n_list1 + n_list2
		if len(n_list) > 0:
			n_dic = collections.Counter(n_list)
			n_info = [ f'{k}({v})' for (k,v) in n_dic.items() ]
		else:
			n_info = None
			
		self.result = t_info, n_info
		
		
	def get_fetchrange(self):
		def get_start12_end12(sv_type, ter1, ter2, pos1, pos2, dist, fors = self.fors, bacs = self.bacs):
			if sv_type == 'INV':
				ThreeToThree = (ter1=='3' and ter2=='3')
				
				pad_start1 = bacs if ThreeToThree else fors
				pad_end2 = fors if ThreeToThree else bacs
				
				if dist >= bacs and dist < fors:
					pad_end1   = dist - bacs if ThreeToThree else bacs
					pad_start2 = bacs if ThreeToThree else dist - bacs
				elif dist < bacs:
					pad_end1 = pad_start2 = dist/2
				else:
					pad_end1   = fors if ThreeToThree else bacs
					pad_start2 = bacs if ThreeToThree else fors
					
			elif sv_type == 'DEL':
				pad_start1 = pad_end2 = bacs
				if dist < fors:
					pad_end1 = pad_start2 = dist/2
				else:
					pad_end1 = pad_start2 = fors
					
			elif sv_type == 'DUP':
				pad_start1 = pad_end2 = fors
				if dist < bacs:
					pad_end1 = pad_start2 = dist/2
				else:
					pad_end1 = pad_start2 = bacs
				
			elif sv_type == 'TRA':
				pad_start1, pad_end1 = (fors, bacs) if ter1=='5' else (bacs, fors)
				pad_start2, pad_end2 = (fors, bacs) if ter2=='5' else (bacs, fors)
					
			start1 = pos1 - pad_start1
			end1   = pos1 + pad_end1
			start2 = pos2 - pad_start2
			end2   = pos2 + pad_end2
					
			return start1, end1, start2, end2
		
		def get_target_start12_end12(sv_type, ter1, ter2, pos1, pos2, dist, start1, end1, start2, end2, fors = self.fors, bacs = self.bacs):
			if sv_type == 'DEL' or sv_type=='DUP':
				target_start1 = start1
				target_end1   = end1
				target_start2 = start2
				target_end2   = end2
				
			elif sv_type=='INV' or sv_type=='TRA':
				target_start1, target_end1 = ((pos1 - fors), (pos1 + bacs)) if ter1=='5' else ((pos1 - bacs), (pos1 + fors))
				target_start2, target_end2 = ((pos2 - fors), (pos2 + bacs)) if ter2=='5' else ((pos2 - bacs), (pos2 + fors))
					
			return target_start1, target_end1, target_start2, target_end2
		
		# MAIN
		start1, end1, start2, end2 = \
			get_start12_end12(self.sv_type, self.ter1, self.ter2, self.pos1, self.pos2, self.dist)
		target_start1, target_end1, target_start2, target_end2 = \
			get_target_start12_end12(self.sv_type, self.ter1, self.ter2, self.pos1, self.pos2, self.dist, start1, end1, start2, end2)
		
		result = dict()
		result['chr1'] = (self.chr1, start1 - 1, end1)
		result['chr2'] = (self.chr2, start2 - 1, end2)
		result['allresults'] = start1, end1, start2, end2, target_start1, target_end1, target_start2, target_end2
		
		self.fetchrange_dict = result
		
		
	def find_SA_reads(self, chr1, start1, end1, chr2, target_start2, target_end2, rplist, chr_list):
		saINFO=[]
		start1=max(start1,1)
		end1=max(end1,1)
		fetchrange = chr1, start1-1, end1
		
		NR = 0
		for rp in rplist:
			if not self.check_within_fetchrange(rp.read, fetchrange):
				continue
				
			if self.readlimit:
				NR += 1
				if NR > self.r_limit:
					break
		
			if \
			(not rp.read.has_tag('SA')) or \
			rp.read.is_secondary or \
			rp.read.is_supplementary:
				continue

			if len(rp.SA_list) == 0: # added on 210128
				continue
				
			for SAdic in rp.SA_list:
				info_ori = ''
				if \
				SAdic['chr'] == chr2 and \
				SAdic['pos'] >= target_start2 and \
				SAdic['pos'] <= target_end2: #check
					
					tmp = readpair_to_SV(
						rp.read.reference_name, rp.read.reference_start + 1, rp.read.is_reverse, rp.M_range, rp.M_range_len, rp.read_size, rp.chr_n,
						SAdic['chr'], SAdic['pos'], SAdic['is_reverse'], SAdic['M_range'], SAdic['M_range_len'], SAdic['read_size'], SAdic['chr_n']
					)
					
					if tmp == None:
						continue
					else:
						MHLEN, bp1, bp2, terminal1, terminal2, rearr, info_ori = tmp
						
			if info_ori=='rs':
				saINFO.append(
					rp.read.reference_name+':'+str(bp1)+';'+\
					SAdic['chr']+':'+str(bp2)+';'+\
					str(MHLEN)+';'+\
					rearr+';'+\
					terminal1+'to'+terminal2
				)
			elif info_ori=='sr':
				saINFO.append(
					SAdic['chr']+':'+str(bp2)+';'+\
					rp.read.reference_name+':'+str(bp1)+';'+\
					str(MHLEN)+';'+\
					rearr+';'+\
					terminal2+'to'+terminal1
				)
		return(saINFO)
			
			
			
# step 03
class calc_final_count(worker):
	# args: chr1, pos1, ter1, chr2, pos2, ter2, rplist_dict, chr_list, r_limit
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.fors = 700 ; self.bacs = 5 ; self.iscut = 700 ; self.sc_co = 5
		
		
	def main(self):
		
		self.get_fetchrange()
		
		sa_seq_list=[]
		t1_list = self.find_discordant_reads_main(
			self.chr1, self.pos1, self.ter1, self.chr2, self.pos2, self.ter2, self.rplist_dict['tumor']['chr1'], sa_seq_list, self.fetchrange_dict
		)
		sa_seq_list=t1_list[6]
		n1_list = self.find_discordant_reads_main(
			self.chr1, self.pos1, self.ter1, self.chr2, self.pos2, self.ter2, self.rplist_dict['normal']['chr1'], sa_seq_list, self.fetchrange_dict
		)
		
		sa_seq_list=[]
		t2_list = self.find_discordant_reads_main(
			self.chr2, self.pos2, self.ter2, self.chr1, self.pos1, self.ter1, self.rplist_dict['tumor']['chr2'], sa_seq_list, self.fetchrange_dict
		)
		sa_seq_list=t2_list[6]
		n2_list = self.find_discordant_reads_main(
			self.chr2, self.pos2, self.ter2, self.chr1, self.pos1, self.ter1, self.rplist_dict['normal']['chr2'], sa_seq_list, self.fetchrange_dict
		)
		
		self.result_step6_chr1 = t1_list[9:11]
		self.result_step6_chr2 = t2_list[9:11]
		
		self.main0(t1_list, n1_list, t2_list, n2_list)
		
		
	def get_fetchrange(self):
		result = dict()
		
		pos1_start, pos1_end, pos2_start, pos2_end = \
			find_discordant_reads_poscalc(self.chr1, self.pos1, self.ter1, self.chr2, self.pos2, self.ter2, self.fors, self.bacs)
		result['chr1'] = self.chr1, pos1_start - 1, pos1_end
		result['chr1_allresults'] = pos1_start, pos1_end, pos2_start, pos2_end
		
		pos1_start, pos1_end, pos2_start, pos2_end = \
			find_discordant_reads_poscalc(self.chr2, self.pos2, self.ter2, self.chr1, self.pos1, self.ter1, self.fors, self.bacs)
		result['chr2'] = self.chr2, pos1_start - 1, pos1_end
		result['chr2_allresults'] = pos1_start, pos1_end, pos2_start, pos2_end
		
		self.fetchrange_dict = result
			
			
	def find_mate_from_SA(self, rp, chr_list):
		
		def find_interCigar_BP(
			chr1, pos1, isrev1, M_range1, M_range_len1, read_size1, chr_n1, 
			chr2, pos2, isrev2, M_range2, M_range_len2, read_size2, chr_n2
		):
			tmp = readpair_to_SV(
				chr1, pos1, isrev1, M_range1, M_range_len1, read_size1, chr_n1, 
				chr2, pos2, isrev2, M_range2, M_range_len2, read_size2, chr_n2
			)
			if tmp == None:
				return 'overlap'
			else:
				MHLEN, bp1, bp2, terminal1, terminal2, rearr, ori = tmp
				info = f'{chr1}:{bp1};{chr2}:{bp2};{MHLEN};{rearr};{terminal1}to{terminal2};{M_range_len1};{M_range_len2}'
				rvs_info = f'{chr2}:{bp2};{chr1}:{bp1};{MHLEN};{rearr};{terminal2}to{terminal1};{M_range_len2};{M_range_len1}'
						
				if chr_n1 < chr_n2: 
					return(info)
				elif chr_n1 > chr_n2:
					return(rvs_info)
				elif chr_n1 == chr_n2:
					if bp1 <= bp2:
						return(info)
					elif bp1 > bp2:
						return(rvs_info)
					
					
		newBP_list=[];neoBP_list=[]
		for SAdic in rp.SA_list:
			res = find_interCigar_BP(
				rp.read.reference_name, rp.read.reference_start+1, rp.read.is_reverse, rp.M_range, rp.M_range_len, rp.read_size, rp.chr_n,
				SAdic['chr'], SAdic['pos'], SAdic['is_reverse'], SAdic['M_range'], SAdic['M_range_len'], SAdic['read_size'], SAdic['chr_n']
			)
			if res != 'overlap':
				newBP_list.append(res)
		if len(rp.SA_list)>1:
			for SAdic1, SAdic2 in itertools.combinations(rp.SA_list,2):
				res = find_interCigar_BP(
					SAdic1['chr'], SAdic1['pos'], SAdic1['is_reverse'], SAdic1['M_range'], SAdic1['M_range_len'], SAdic1['read_size'], SAdic1['chr_n'],
					SAdic2['chr'], SAdic2['pos'], SAdic2['is_reverse'], SAdic2['M_range'], SAdic2['M_range_len'], SAdic2['read_size'], SAdic2['chr_n']
				)
				if res != 'overlap':
					neoBP_list.append(res)
		return(newBP_list, neoBP_list)
		
		
	def find_discordant_reads_main(self, chr1, pos1, ter1, chr2, pos2, ter2, rplist, sa_seq_list, fetchrange_dict):
		
		def mate_list_summary(mate_list):
			summary_dic={}
			for mate in mate_list:
				mate_indi=mate.split(';')
				m1=int(mate_indi[5])
				m2=int(mate_indi[6])
				info=';'.join(mate_indi[0:5])
				if info not in summary_dic.keys():
					summary_dic[info]={}
					summary_dic[info]['num']=0
					summary_dic[info]['match1']=[]
					summary_dic[info]['match2']=[]
				summary_dic[info]['num']+=1
				summary_dic[info]['match1'].append(m1)
				summary_dic[info]['match2'].append(m2)
			final_list=[]
			for info in summary_dic.keys():
				m1max=max(summary_dic[info]['match1'])
				m2max=max(summary_dic[info]['match2'])
				freq=summary_dic[info]['num']
				final_list.append(info+';'+str(m1max)+';'+str(m2max)+'('+str(freq)+')')
			return (','.join(final_list))
		
		pair_true_list=[];sp_true_list=[];sa_true_list=[]
		pair_ref_list=[]; jx_ref_list=[]

		new_mate_list=[];neo_mate_list=[]

		sa_seq_list_internal = []
		true_mapq_list=[]; true_pos_list=[]
		
		if chr1 == self.chr1:
			fetchrange = fetchrange_dict['chr1']
			pos1_start, pos1_end, pos2_start, pos2_end = fetchrange_dict['chr1_allresults']
		elif chr1 == self.chr2:
			fetchrange = fetchrange_dict['chr2']
			pos1_start, pos1_end, pos2_start, pos2_end = fetchrange_dict['chr2_allresults']
		
		NR = 0
		for rp in rplist:
			if not self.check_within_fetchrange(rp.read, fetchrange):
				continue
				
			if self.readlimit:
				NR += 1
				if NR > self.r_limit:
					break
				
			if \
			(not rp.read.is_paired) or \
			rp.read.mate_is_unmapped:
				continue
			
			if \
			(ter1=='3' and (rp.read.cigartuples[-1][0]==4 or rp.read.cigartuples[-1][0]==5)) or \
			(ter1=='5' and (rp.read.cigartuples[0][0] ==4 or rp.read.cigartuples[0][0] ==5)):
				flag = True
			else:
				flag = False
				
			for SAdic in rp.SA_list:
				if readfilter_1(SAdic, chr2, pos2):
					true_mapq_list.append(rp.read.mapping_quality)
					
					if rp.read.cigartuples[0][0]==4 or rp.read.cigartuples[0][0]==5:
						true_pos_list.append(rp.read.reference_start+1-rp.read.cigartuples[0][1])
					else:
						true_pos_list.append(rp.read.reference_start+1)
					
					if flag:
						sa_true_list.append(rp.read.query_name)
						sp_true_list.append(rp.read.query_name)

					if ter1=='3' and rp.read.cigartuples[-1][0]==4:
						#sc_seq=rp.read.query_sequence[rp.read.cigartuples[-1][1]*(-1) : rp.read.cigartuples[-1][1]*(-1)+self.sc_co]
						sa_seq_list_internal.append(rp.sc_seq)
						
						if flag:
							#sc_seq=rp.read.query_sequence[rp.read.cigartuples[-1][1]*(-1): rp.read.cigartuples[-1][1]*(-1)+self.sc_co]
							sa_seq_list.append(rp.sc_seq)
							sa_res=self.find_mate_from_SA(rp, self.chr_list)
							new_mate_list += sa_res[0]
							neo_mate_list += sa_res[1]
					elif ter1=='5' and rp.read.cigartuples[0][0]==4:
						#sc_seq=rp.read.query_sequence[rp.read.cigartuples[0][1]-1-self.sc_co+1 : rp.read.cigartuples[0][1]-1+1]
						sa_seq_list_internal.append(rp.sc_seq)
						
						if flag:
							#sc_seq=rp.read.query_sequence[rp.read.cigartuples[0][1]-1-self.sc_co+1 : rp.read.cigartuples[0][1]-1+1]
							sa_seq_list.append(rp.sc_seq)
							sa_res=self.find_mate_from_SA(rp, self.chr_list)
							new_mate_list += sa_res[0]
							neo_mate_list += sa_res[1]
						
		sa_seq_list=list(set(sa_seq_list))
		sa_seq_list_internal=list(set(sa_seq_list_internal))
		
		
		NR = 0
		for rp in rplist:
			if not self.check_within_fetchrange(rp.read, fetchrange):
				continue
					
			if self.readlimit:
				NR += 1
				if NR > self.r_limit:
					break
				
			if \
			(not rp.read.is_paired) or \
			rp.read.mate_is_unmapped:
				continue
				
			pair_ref_mode='off';#jx_ref_mode='off'
			if ter1=='3':
				if \
				(not rp.read.is_reverse) and \
				rp.read.mate_is_reverse and \
				rp.read.next_reference_name == chr1 and \
				rp.read.reference_start +1 < pos1 and \
				rp.read.reference_start +1 + rp.read.template_length -1 > pos1 and \
				rp.read.template_length >= 0 and \
				rp.read.template_length < self.iscut: 
					pair_ref_list.append(rp.read.query_name)
					pair_ref_mode='on'
					
				if \
				rp.read.reference_start + 1 <= pos1 and \
				rp.read.reference_start + 1 + rp.read.reference_length - 1 > pos1 and \
				rp.read.next_reference_name == chr1:
					jx_ref_list.append(rp.read.query_name)
					#jx_ref_mode='on'
					
				if \
				pair_ref_mode == 'off' and \
				(not rp.read.is_reverse) and \
				rp.read.next_reference_name == chr2 and \
				rp.read.next_reference_start +1 >= pos2_start and \
				rp.read.next_reference_start +1 < pos2_end:
					
					if \
					(ter2=='3' and (not rp.read.mate_is_reverse)) or \
					(ter2=='5' and rp.read.mate_is_reverse): 
						
						true_mapq_list.append(rp.read.mapping_quality)
						if rp.read.cigartuples[0][0]==4 or rp.read.cigartuples[0][0]==5:
							true_pos_list.append(rp.read.reference_start+1-rp.read.cigartuples[0][1])
						else:
							true_pos_list.append(rp.read.reference_start+1)
						
						pair_true_list.append(rp.read.query_name)
						if \
						rp.read.has_tag('SA') and \
						rp.read.query_name not in sa_true_list and \
						rp.read.query_name not in pair_ref_list:
							sa_res=self.find_mate_from_SA(rp, self.chr_list)
							new_mate_list += sa_res[0]
							neo_mate_list += sa_res[1]
							
				if \
				pos1 - (rp.read.reference_start + 1) + 1 == rp.read.reference_length and \
				(rp.read.cigartuples[-1][0]==4 and rp.read.cigartuples[-1][1] >= self.sc_co):
					#sc_seq=rp.read.query_sequence[rp.read.cigartuples[-1][1]*(-1): rp.read.cigartuples[-1][1]*(-1)+self.sc_co]
							
					if rp.sc_seq in sa_seq_list:
						sp_true_list.append(rp.read.query_name)
						
					if rp.sc_seq in sa_seq_list_internal:
						true_mapq_list.append(rp.read.mapping_quality)
						if rp.read.cigartuples[0][0]==4 or rp.read.cigartuples[0][0]==5:
							true_pos_list.append(rp.read.reference_start+1-rp.read.cigartuples[0][1])
						else:
							true_pos_list.append(rp.read.reference_start+1)
								
								
			elif ter1=='5':
				if \
				rp.read.is_reverse and \
				(not rp.read.mate_is_reverse) and \
				rp.read.next_reference_name==chr1 and \
				rp.read.reference_start +1 + rp.read.reference_length -1 >= pos1 and \
				rp.read.reference_start + 1 + rp.read.reference_length -1 + rp.read.template_length + 1 < pos1 and \
				rp.read.template_length < 0 and \
				rp.read.template_length*(-1) < self.iscut:  # in this situation read.template_length is negative value
					pair_ref_list.append(rp.read.query_name)
					pair_ref_mode='on'
					
				if \
				rp.read.reference_start + 1 < pos1 and \
				rp.read.reference_start + 1 + rp.read.reference_length - 1 >= pos1: 
					jx_ref_list.append(rp.read.query_name)
					#jx_ref_mode='on'
					
				if \
				pair_ref_mode=='off' and \
				rp.read.is_reverse and \
				rp.read.next_reference_name == chr2 and \
				rp.read.next_reference_start +1 >= pos2_start and \
				rp.read.next_reference_start +1 < pos2_end:
					if (ter2=='3' and (not rp.read.mate_is_reverse)) or (ter2=='5' and rp.read.mate_is_reverse):
						
						true_mapq_list.append(rp.read.mapping_quality)
						if rp.read.cigartuples[0][0]==4 or rp.read.cigartuples[0][0]==5:
							true_pos_list.append(rp.read.reference_start+1-rp.read.cigartuples[0][1])
						else:
							true_pos_list.append(rp.read.reference_start+1)
						
						pair_true_list.append(rp.read.query_name)
						if \
						rp.read.has_tag('SA') and \
						rp.read.query_name not in sa_true_list and \
						rp.read.query_name not in pair_ref_list:
							sa_res=self.find_mate_from_SA(rp, self.chr_list)
							new_mate_list += sa_res[0]
							neo_mate_list += sa_res[1]
							
				if \
				rp.read.reference_start + 1 == pos1 and \
				(rp.read.cigartuples[0][0] == 4 and rp.read.cigartuples[0][1] >= self.sc_co):
					#sc_seq = rp.read.query_sequence[rp.read.cigartuples[0][1]-1-self.sc_co+1 : rp.read.cigartuples[0][1]-1+1]
					
					if rp.sc_seq in sa_seq_list:
						sp_true_list.append(rp.read.query_name)
						
					if rp.sc_seq in sa_seq_list_internal:
						true_mapq_list.append(rp.read.mapping_quality)
						if rp.read.cigartuples[0][0]==4 or rp.read.cigartuples[0][0]==5:
							true_pos_list.append(rp.read.reference_start+1-rp.read.cigartuples[0][1])
						else:
							true_pos_list.append(rp.read.reference_start+1)
								
		sa_true_list=list(set(sa_true_list))
		pair_ref_list=list(set(pair_ref_list))
		jx_ref_list=list(set(jx_ref_list) & set(pair_ref_list))
		all_ref_list=list(set(pair_ref_list+jx_ref_list)-set(sa_true_list))
		pair_true_list=list(set(pair_true_list)-set(all_ref_list))
		sp_true_list=list(set(sp_true_list))
		all_true_list=list(set(pair_true_list+sp_true_list+sa_true_list))
		
		if len(new_mate_list)==0:
			new_mate='NA'
		else:
			new_mate=mate_list_summary(new_mate_list)
			
		if len(neo_mate_list)==0:
			neo_mate='NA'
		else:
			neo_mate=mate_list_summary(neo_mate_list)
			
		return([
			pair_true_list, 
			sp_true_list, 
			sa_true_list, 
			pair_ref_list, 
			jx_ref_list, 
			all_ref_list, 
			sa_seq_list, 
			new_mate, 
			neo_mate,
			true_mapq_list,
			true_pos_list,
		])
		
	
	def main0(self, t1_list, n1_list, t2_list, n2_list):
		a1=0; as1=0; asa1=0;r1=0;rj1=0; r2=0; rj2=0; na1=0; nsa1=0
		#normal_split1='off';normal_split2='off'
		
		t1_pair_list=t1_list[0]
		t1_sp_list=t1_list[1]
		t1_sa_list=t1_list[2]
		t1_rj_list=t1_list[4]
		t1_rt_list=t1_list[5]
		new_mate1=t1_list[7]
		neo_mate1=t1_list[8]
		
		n1_pair_list=n1_list[0]
		n1_sp_list=n1_list[1]
		n1_sa_list=n1_list[2]
		n1_rt_list=n1_list[5]
		
		t2_pair_list=t2_list[0]
		t2_sp_list=t2_list[1]
		t2_sa_list=t2_list[2]
		t2_rj_list=t2_list[4]
		t2_rt_list=t2_list[5]
		new_mate2=t2_list[7]
		neo_mate2=t2_list[8]
		
		n2_pair_list=n2_list[0]
		n2_sp_list=n2_list[1]
		n2_sa_list=n2_list[2]
		n2_rt_list=n2_list[5]
		
		t1_total_list=list(set(t1_pair_list+t1_sp_list+t1_sa_list))
		t1_sp_list=list(set(t1_sp_list+t1_sa_list))
		t2_total_list=list(set(t2_pair_list+t2_sp_list+t2_sa_list))
		t2_sp_list=list(set(t2_sp_list+t2_sa_list))
		n1_sp_list=list(set(n1_sp_list+n1_sa_list))
		n2_sp_list=list(set(n2_sp_list+n2_sa_list))
		n_split_n=len(list(set(n1_sp_list+n2_sp_list)))

		t_tot_n=len(list(set(t1_total_list+t2_total_list)))
		t_split_n=len(list(set(t1_sp_list+t2_sp_list)))
		t_sa_n=len(list(set(t1_sa_list+t2_sa_list)))
		t1_reftot_n=len(t1_rt_list)
		t1_refjx_n=len(t1_rj_list)
		t2_reftot_n=len(t2_rt_list)
		t2_refjx_n=len(t2_rj_list)
		n1_reftot_n=len(n1_rt_list)
		n2_reftot_n=len(n2_rt_list)
		n_tot_n=len(list(set(n1_pair_list+n1_sa_list+n2_pair_list+n2_sa_list)))
		n_sa_n=len(list(set(n1_sa_list+n2_sa_list)))
		
		self.result = [
			t_tot_n, t_split_n, t_sa_n, t1_reftot_n, t1_refjx_n, t2_reftot_n, t2_refjx_n, n_tot_n, n_sa_n, 
			new_mate1, neo_mate1, new_mate2, neo_mate2, n1_reftot_n, n2_reftot_n, n_split_n
		]
	
		
class count_frag_num(worker):
	# args: chr1, pos1, rplist, r_limit
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		
	def get_fetchrange(self):
		self.fetchrange = (self.chr1, self.pos1 - 1, self.pos1)
	
	def main(self):
		self.get_fetchrange()
		
		total_frag_list=[]
		
		NR = 0
		for rp in self.rplist:
			if not self.check_within_fetchrange(rp.read, self.fetchrange):
				continue
				
			if self.readlimit:
				NR += 1
				if NR > self.r_limit:
					break
				
			if \
			(not rp.read.is_paired) or \
			rp.read.mate_is_unmapped:
				continue
			total_frag_list.append(rp.read.query_name)
			
		total_frag_list=list(set(total_frag_list))
		self.result = len(total_frag_list)
	
	
# step 4	
class find_pnsc(worker):
	# args: chr1, pos1, ter1, chr2, pos2, ter2, tumor_rplist, normal_rplist, r_limit
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.fors = 500 ; self.bacs = 5 ; self.sc_co = 5 ; self.n_sr = 10
		
	def get_fetchrange(self):
		result = dict()
		pos1_start, pos1_end, pos2_start, pos2_end = \
			find_discordant_reads_poscalc(self.chr1, self.pos1, self.ter1, self.chr2, self.pos2, self.ter2, self.fors, self.bacs)
		result['tumor'] = ( self.chr1, pos1_start - 1, pos1_end )
		result['normal'] = ( self.chr1, max(1,self.pos1-1-self.n_sr), self.pos1+self.n_sr )
		result['allresults'] = pos1_start, pos1_end, pos2_start, pos2_end
		self.fetchrange_dict = result
	
	def main(self):
		self.get_fetchrange()
		sa_seq_list=[];sp_true_list=[]
		pos1_start, pos1_end, pos2_start, pos2_end = self.fetchrange_dict['allresults']
		
		NR = 0
		for rp in self.tumor_rplist:
			if not self.check_within_fetchrange(rp.read, self.fetchrange_dict['tumor']):
				continue
				
			if self.readlimit:
				NR += 1
				if NR > self.r_limit:
					break
		
			if \
			(not rp.read.is_paired) or \
			rp.read.mate_is_unmapped or \
			rp.read.is_secondary or \
			rp.read.is_supplementary:
				continue
				
			for SAdic in rp.SA_list:
				if readfilter_1(SAdic, self.chr2, self.pos2):
					if self.ter1=='3' and rp.read.cigartuples[-1][0]==4 and rp.read.cigartuples[-1][1] >= self.sc_co:
						sa_seq_list.append(rp.sc_seq)
					elif self.ter1=='5' and rp.read.cigartuples[0][0]==4 and rp.read.cigartuples[0][1] >= self.sc_co:
						sa_seq_list.append(rp.sc_seq)
							
			if self.ter1=='3':
				if \
				(not rp.read.is_reverse) and \
				rp.read.next_reference_name == self.chr2 and \
				rp.read.next_reference_start +1 >= pos2_start and \
				rp.read.next_reference_start +1 < pos2_end:
					if \
					(self.ter2=='3' and (not rp.read.mate_is_reverse)) or \
					(self.ter2=='5' and rp.read.mate_is_reverse): 
						if rp.read.cigartuples[-1][0]==4 and rp.read.cigartuples[-1][1] >= self.sc_co:
							sa_seq_list.append(rp.sc_seq)
			elif self.ter1=='5':
				if \
				rp.read.is_reverse and \
				rp.read.next_reference_name == self.chr2 and \
				rp.read.next_reference_start +1 >= pos2_start and \
				rp.read.next_reference_start +1 < pos2_end:
					if \
					(self.ter2=='3' and (not rp.read.mate_is_reverse)) or \
					(self.ter2=='5' and rp.read.mate_is_reverse):
						if rp.read.cigartuples[0][0]==4 and rp.read.cigartuples[0][1] >= self.sc_co:
							sa_seq_list.append(rp.sc_seq)
							
		sa_seq_list=list(set(sa_seq_list))
		
		
		NR = 0
		for rp in self.normal_rplist:
			if not self.check_within_fetchrange(rp.read, self.fetchrange_dict['normal']):
				continue
				
			if self.readlimit:
				NR += 1
				if NR > self.r_limit:
					break
				
			if \
			(not rp.read.is_paired) or \
			rp.read.mate_is_unmapped or \
			rp.read.is_secondary or \
			rp.read.is_supplementary:
				continue
				
			if self.ter1=='3':
				if len(sa_seq_list) > 0:
					if rp.read.cigartuples[-1][0]==4 and rp.read.cigartuples[-1][1] >= self.sc_co:
						if rp.sc_seq in sa_seq_list:
							sp_true_list.append(rp.read.query_name)
			elif self.ter1=='5':
				if len(sa_seq_list) > 0:
					if rp.read.cigartuples[0][0] == 4 and rp.read.cigartuples[0][1] >= self.sc_co:
						if rp.sc_seq in sa_seq_list:
							sp_true_list.append(rp.read.query_name)
							
		self.result = len(list(set(sp_true_list)))
		
		
# step 5
class amount_discordant(worker):
	# args: chr1, pos1, rplist, r_limit
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.ser = 300 ; self.iscut = 1000 ; self.mate_bin = 100000
		
	def get_fetchrange(self):
		pos1_start = max(self.pos1 - self.ser, 1)
		pos1_end = self.pos1 + self.ser
		self.fetchrange = self.chr1, pos1_start-1, pos1_end
		
	def normal_fragment_filter(self, rp):
		if \
		(
			rp.read.next_reference_name == self.chr1 and \
			(not rp.read.is_reverse) and \
			rp.read.template_length > 0 and \
			rp.read.template_length < self.iscut \
		) or \
		(
			rp.read.next_reference_name == self.chr1 and \
			rp.read.is_reverse and \
			rp.read.template_length < 0 and \
			rp.read.template_length*(-1) < self.iscut
		):
			return True
		else:
			return False
	
	def main(self):
		disc_mate_dic={};normal_frag_list=[];total_frag_list=[];disc_frag_list=[];tra_frag_list=[];clip_frag_list=[]
		
		self.get_fetchrange()
		
		NR = 0
		for rp in self.rplist:
			if not self.check_within_fetchrange(rp.read, self.fetchrange):
				continue
				
			if self.readlimit:
				NR += 1
				if NR > self.r_limit:
					break
				
			if \
			(not rp.read.is_paired) or \
			rp.read.mate_is_unmapped or \
			rp.read.mapping_quality < 1:
				continue
				
				
			total_frag_list.append(rp.read.query_name)
			
			if sum( rp.read.get_cigar_stats()[1][4:6] ) > 0:
				clip_frag_list.append(rp.read.query_name)
				
			if self.normal_fragment_filter(rp):
				normal_frag_list.append(rp.read.query_name)
			else:
				disc_frag_list.append(rp.read.query_name)
				
				if rp.read.next_reference_name != self.chr1:
					tra_frag_list.append(rp.read.query_name)
					
				if rp.read.next_reference_name not in disc_mate_dic:
					disc_mate_dic[rp.read.next_reference_name] = dict()
				binn=(rp.read.next_reference_start+1)/self.mate_bin
				if binn not in disc_mate_dic[rp.read.next_reference_name]:
					disc_mate_dic[rp.read.next_reference_name][binn] = list()
				disc_mate_dic[rp.read.next_reference_name][binn].append(rp.read.query_name)
		
		total_fn=len(list(set(total_frag_list)))
		normalp_fn=len(list(set(normal_frag_list)))
		disc_fn=len(list(set(disc_frag_list)))
		tra_fn=len(list(set(tra_frag_list)))
		clip_fn=len(list(set(clip_frag_list)))
		disc_chr_n=len(disc_mate_dic)
		
		disc_bin_n=0;disc_bin_n2=0
		disc_chr_n2_list=[]
		for chrom in disc_mate_dic:
			disc_bin_n += len(disc_mate_dic[chrom])
			if len(disc_mate_dic[chrom]) >= 2 :
				disc_chr_n2_list.append(chrom)
			for eachbin in disc_mate_dic[chrom]:
				if len(disc_mate_dic[chrom][eachbin]) >=2:
					disc_bin_n2 +=1
					disc_chr_n2_list.append(chrom)
		disc_chr_n2=len(list(set(disc_chr_n2_list)))
		
		info_list=[
			str(self.ser), str(self.iscut), str(round(self.mate_bin/float(1000),1))+'k', str(total_fn), str(normalp_fn), str(disc_fn), str(tra_fn),
			str(disc_chr_n), str(disc_bin_n), str(disc_chr_n2), str(disc_bin_n2)
		]
		self.result = ';'.join(info_list) + '\t' + str(clip_fn) + '\t' + str(disc_bin_n)
		
def set_constants(refver):
	# set satellite regions and primary contigs
	satellite_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
	if refver == '19':
		satellite_path = satellite_path + '/data/excluded_region/satellite_hg19.bed.gz'
		primary_contigs = [ str(x) for x in list(range(1, 23)) + ['X', 'Y'] ]
	elif refver == '38':
		satellite_path = satellite_path + '/data/excluded_region/satellite_hg38.bed.gz'
		primary_contigs = [ 'chr'+str(x) for x in list(range(1, 23)) + ['X', 'Y'] ]

	pr_satellite = pr.read_bed(satellite_path)

	# others
	r_limit = 50000 # The number of reads that this script searches for analysis
	r_limit_dp = 1000 # In the initial fetch step, average depth higher than this will be trimmed
	terinfo_list = ['3to3', '3to5', '5to3', '5to5']

	return satellite_path, pr_satellite, primary_contigs, r_limit, r_limit_dp, terinfo_list


def printErr(*args, **kwargs):
	print(*args, file = sys.stderr, flush = True, **kwargs)

	
# reference chr size (for step 7)
def get_chr_size(fai_path):
	chr_size = dict()
	chr_list = list()
	with open(fai_path, 'r') as ref_file:
		for line in ref_file:
			linesp = line.split('\t')
			chrom = linesp[0] ; length = int(linesp[1])
			chr_size[chrom] = length
			chr_list.append(chrom)

	return chr_size, chr_list


def get_header(firstline):
	hdr = firstline.strip().split('\t')
	hdr.append('tBPinfo')
	hdr.append('nBPinfo')
	hdr.append('re_chr1')
	hdr.append('re_pos1')
	hdr.append('re_chr2')
	hdr.append('re_pos2')
	hdr.append('MH')
	hdr.append('terminal')
	hdr.append('SVtype')
	hdr.append('tSA')
	hdr.append('nSA')
	hdr.append(';'.join([
		'Tumor_Ref1', 
		'Ref2', 
		'AllDiscordantFragments', 
		'SplitFragments', 
		'SATagFragments', 
		'Vaf1', 
		'Vaf2'
	]))
	hdr.append(';'.join([
		'PairNormal_Ref1', 
		'Ref2', 
		'AllDiscordantFragments', 
		'SplitFragments', 
		'SATagFragments', 
		'FragCount1', 
		'FragCount2'
	]))
	hdr.append('new_mate1')
	hdr.append('neo_mate1')
	hdr.append('new_mate2')
	hdr.append('neo_mate2')
	hdr.append('PairNormalSameClip')
	hdr.append(';'.join([
		'T_BP1_SearchRange', 
		'NormalDistance', 
		'MateBin', 
		'TotalFragN', 
		'NormalPairFragN', 
		'DiscordantFragN', 
		'TRAFragN', 
		'DiscordantChrN', 
		'DiscordantBinN', 
		'DiscordantChrN2', 
		'DiscordantBinN2'
	]))
	hdr.append('T_BP1_clip_readN')
	hdr.append('T_BP1_other_discordant_cluster')
	hdr.append(';'.join([
		'T_BP2_SearchRange', 
		'NormalDistance', 
		'MateBin', 
		'TotalFragN', 
		'NormalPairFragN', 
		'DiscordantFragN', 
		'TRAFragN', 
		'DiscordantChrN', 
		'DiscordantBinN', 
		'DiscordantChrN2', 
		'DiscordantBinN2'
	]))
	hdr.append('T_BP2_clip_readN')
	hdr.append('T_BP2_other_discordant_cluster')
	hdr.append(';'.join([
		'N_BP1_SearchRange',
		'NormalDistance',
		'MateBin',
		'TotalFragN',
		'NormalPairFragN',
		'DiscordantFragN',
		'TRAFragN',
		'DiscordantChrN',
		'DiscordantBinN',
		'DiscordantChrN2',
		'DiscordantBinN2'
	]))
	hdr.append('N_BP1_clip_readN')
	hdr.append('N_BP1_other_discordant_cluster')
	hdr.append(';'.join([
		'N_BP2_SearchRange',
		'NormalDistance',
		'MateBin',
		'TotalFragN',
		'NormalPairFragN',
		'DiscordantFragN',
		'TRAFragN',
		'DiscordantChrN',
		'DiscordantBinN',
		'DiscordantChrN2',
		'DiscordantBinN2'
	]))
	hdr.append('N_BP2_clip_readN')
	hdr.append('N_BP2_other_discordant_cluster')
	hdr.append('MAPQ1_min;med;max')
	hdr.append('MAPQ2_min;med;max')
	hdr.append('POS1_min;med;max')
	hdr.append('POS2_min;med;max')
	hdr.append('depth_ratio_change_bp1')
	hdr.append('depth_ratio_change_bp2')
	
	return hdr


def step1(chr1, pos1, ter1, chr2, pos2, ter2, sv_type, dist, chr_list, rplist_dict, r_limit):
	wrapper_step1 = step1_wrapper(chr1=chr1, pos1=pos1, ter1=ter1, chr2=chr2, pos2=pos2, ter2=ter2, 
									sv_type=sv_type, dist=dist, chr_list=chr_list, rplist_dict=rplist_dict,
									r_limit=r_limit)
	wrapper_step1.main()
	return wrapper_step1.result

def step2(chr1, pos1, chr2, pos2, terinfo, sv_type, t_info, n_info, line):
	def BP_adjustment(info, terinfo, sv_type, chr1, pos1, chr2, pos2, line):
		def subr(BP):
			BP = BP.replace('NC_','NC').split(';')
			BP_sp0 = BP[0].split(':')
			BP_sp1 = BP[1].split(':')
			
			bpchr1 = BP_sp0[0].replace('NC','NC_')
			bppos1 = int( BP_sp0[1] )
			bpchr2 = BP_sp1[0].replace('NC','NC_')
			bppos2 = int( BP_sp1[1] )
			bpMH = BP[2] 
			
			if chr1!=bpchr1:
				printErr('chromosome1 error')
				printErr(candidate_list)
				sys.exit(1)
			if chr2!=bpchr2:
				printErr('chromosome2 error')
				printErr(line, end='')
				sys.exit(1)
			
			return bpchr1, bppos1, bpchr2, bppos2, bpMH
		
		
		tSA = 0 ; MH = '.' # initial value
		
		if info == None:
			return chr1, pos1, chr2, pos2, MH, tSA
		else:
			current_tSA=0
			candidate_list=[]
			for BP in info:
				if terinfo in BP and sv_type in BP:
					tSA_count=int(BP.split('(')[1][:-1])
					if tSA_count > current_tSA:
						candidate_list=[]
						candidate_list.append(BP)
						current_tSA=tSA_count
					elif tSA_count == current_tSA:
						candidate_list.append(BP)
					elif tSA_count < current_tSA:
						pass
				
			if len(candidate_list)==1:
				bpchr1, bppos1, bpchr2, bppos2, MH = subr(candidate_list[0])
				pos1 = bppos1
				pos2 = bppos2
				tSA = current_tSA
				
			elif len(candidate_list) >1:
				current_distsum = None
				for BP in candidate_list:
					bpchr1, bppos1, bpchr2, bppos2, MH = subr(BP)
					distsum = abs(pos1-bppos1) + abs(pos2-bppos2)
					if current_distsum == None or distsum < current_distsum:
						final_info=[bpchr1, bppos1, bpchr2, bppos2, MH]
						current_distsum=distsum
						
				chr1, pos1, chr2, pos2, MH = final_info
				tSA = current_tSA
				
			return chr1, pos1, chr2, pos2, MH, tSA
		
		
	def get_nSA(n_info, terinfo, sv_type, pos1, pos2):
		nSA = 0 # initial value
		
		if n_info == None:
			return nSA
		else:
			nSA_candidate=[]
			for BP in n_info:
				if terinfo in BP and sv_type in BP:
					BP = BP.replace('NC_','NC')
					BPsp = BP.split(';')
					bppos1=int((BPsp[0]).split(':')[1])
					bppos2=int((BPsp[1]).split(':')[1])
					nSA_count=int(BP.split('(')[1][:-1])

					dist1=abs(pos1-bppos1)
					dist2=abs(pos2-bppos2)
					if dist1 <= 1 and dist2 <=1:
						nSA_candidate.append(nSA_count)
			if len(nSA_candidate)==0:
				nSA = 0
			else:
				nSA = max(nSA_candidate)
				
			return nSA
		
	# MAIN
	chr1, pos1, chr2, pos2, MH, tSA = BP_adjustment(t_info, terinfo, sv_type, chr1, pos1, chr2, pos2, line)
	nSA = get_nSA(n_info, terinfo, sv_type, pos1, pos2)
	info_list = [chr1, str(pos1), chr2, str(pos2), MH, terinfo, sv_type, str(tSA), str(nSA)]
	
	return chr1, pos1, chr2, pos2, info_list
	
	
def step3(chr1, pos1, ter1, chr2, pos2, ter2, rplist_dict, wrapper_calcFinalCount, r_limit, shortDco = 500):
	if pos2 == '.': 
		print_list = [line.strip()] + ['.']*6
	else:
		adf, \
		sf, \
		saf, \
		ref1, \
		rj1, \
		ref2, \
		rj2, \
		na1, \
		nsa1, \
		new_mate1, \
		neo_mate1, \
		new_mate2, \
		neo_mate2, \
		n_ref1, \
		n_ref2, \
		nsf = wrapper_calcFinalCount.result
		
		wrapper_countFragNum = count_frag_num(chr1=chr1, pos1=pos1, rplist=rplist_dict['normal']['chr1'], r_limit=r_limit)
		wrapper_countFragNum.main()
		pnfc1 = wrapper_countFragNum.result
		
		wrapper_countFragNum = count_frag_num(chr1=chr2, pos1=pos2, rplist=rplist_dict['normal']['chr2'], r_limit=r_limit)
		wrapper_countFragNum.main()
		pnfc2 = wrapper_countFragNum.result

		if chr1 == chr2 and ter1 == '3' and ter2 == '5' and abs(pos2-pos1) < shortDco:
			adf=sf
			ref1=rj1
			ref2=rj2
		elif adf == sf:
			ref1=rj1
			ref2=rj2
			
		vaf1 = 'NA' if (adf+ref1 == 0) else str(round((adf)*100/float(adf+ref1),2))+'%'
		vaf2 = 'NA' if (adf+ref2 == 0) else str(round((adf)*100/float(adf+ref2),2))+'%'
			
		# asr1 and asr2 were not counted in 'else' d/t redundancy with r1, r2
		t_info_list=[str(ref1),str(ref2),str(adf),str(sf),str(saf),vaf1,vaf2]
		n_info_list=[str(n_ref1), str(n_ref2), str(na1), str(nsf), str(nsa1),str(pnfc1), str(pnfc2)]
		print_list=[';'.join(t_info_list), ';'.join(n_info_list), new_mate1, neo_mate1, new_mate2, neo_mate2]
		
	return print_list


def step4(chr1, pos1, ter1, chr2, pos2, ter2, rplist_dict, r_limit):
	
	wrapper_findPnsc = find_pnsc(
		chr1=chr1, pos1=pos1, ter1=ter1, chr2=chr2, pos2=pos2, ter2=ter2,
		tumor_rplist=rplist_dict['tumor']['chr1'], normal_rplist=rplist_dict['normal']['chr1'],
		r_limit=r_limit
	)
	wrapper_findPnsc.main()
	pnsc1 = wrapper_findPnsc.result
	
	wrapper_findPnsc = find_pnsc(
		chr1=chr2, pos1=pos2, ter1=ter2, chr2=chr1, pos2=pos1, ter2=ter1,
		tumor_rplist=rplist_dict['tumor']['chr2'], normal_rplist=rplist_dict['normal']['chr2'],
		r_limit=r_limit
	)
	wrapper_findPnsc.main()
	pnsc2 = wrapper_findPnsc.result
	
	pnsc = pnsc1 + pnsc2
	
	return pnsc


def step5(chr1, pos1, chr2, pos2, rplist_dict, r_limit):
	
	wrapper_amountDiscordant = amount_discordant(chr1=chr1, pos1=pos1, rplist=rplist_dict['tumor']['chr1'], r_limit=r_limit)
	wrapper_amountDiscordant.main()
	tres1 = wrapper_amountDiscordant.result
	
	wrapper_amountDiscordant = amount_discordant(chr1=chr2, pos1=pos2, rplist=rplist_dict['tumor']['chr2'], r_limit=r_limit)
	wrapper_amountDiscordant.main()
	tres2 = wrapper_amountDiscordant.result
	
	wrapper_amountDiscordant = amount_discordant(chr1=chr1, pos1=pos1, rplist=rplist_dict['normal']['chr1'], r_limit=r_limit)
	wrapper_amountDiscordant.main()
	nres1 = wrapper_amountDiscordant.result
	
	wrapper_amountDiscordant = amount_discordant(chr1=chr2, pos1=pos2, rplist=rplist_dict['normal']['chr2'], r_limit=r_limit)
	wrapper_amountDiscordant.main()
	nres2 = wrapper_amountDiscordant.result
	
	return tres1, tres2, nres1, nres2
	
	
def step6(wrapper_calcFinalCount):
	mapq_list1, pos_list1 = wrapper_calcFinalCount.result_step6_chr1
	mapq_list2, pos_list2 = wrapper_calcFinalCount.result_step6_chr2
		
	if len(mapq_list1)==0:
		mq_info1='NA'
	else:
		mq_med1=median(mapq_list1)
		mq_min1=min(mapq_list1)
		mq_max1=max(mapq_list1)
		mq_info1=str(mq_min1)+';'+str(mq_med1)+';'+str(mq_max1)
		
	if len(mapq_list2)==0:
		mq_info2='NA'
	else:
		mq_med2=median(mapq_list2)
		mq_min2=min(mapq_list2)
		mq_max2=max(mapq_list2)
		mq_info2=str(mq_min2)+';'+str(mq_med2)+';'+str(mq_max2)
		
	if len(pos_list1)==0:
		pos_info1='NA'
	else:
		pos_min1=min(pos_list1)
		pos_med1=median(pos_list1)-pos_min1
		pos_max1=max(pos_list1)-pos_min1
		pos_info1='0;'+str(pos_med1)+';'+str(pos_max1)
		
	if len(pos_list2)==0:
		pos_info2='NA'
	else:
		pos_min2=min(pos_list2)
		pos_med2=median(pos_list2)-pos_min2
		pos_max2=max(pos_list2)-pos_min2
		pos_info2='0;'+str(pos_med2)+';'+str(pos_max2)
			
	return mq_info1, mq_info2, pos_info1, pos_info2


def step7(chr1, pos1, ter1, chr2, pos2, ter2, tbam_file, nbam_file, chr_size, searchlen = 500):
	def subr(chrom, pos, bam_file):
		dp_listL=[]
		dp_listR=[]
		startL = max(pos-1-searchlen, 1)
		endL   = max(pos, 2)
		startR = max(pos-1, 1)
		endR   = min(pos+searchlen, chr_size[chrom])

		distL = endL - startL
		distR = endR - startR

		posesL = bam_file.count_coverage(chrom, startL, endL, read_callback='nofilter', quality_threshold = 0)
		posesR = bam_file.count_coverage(chrom, startR, endR, read_callback='nofilter', quality_threshold = 0)
		for idx in range(distL):
			dp_listL.append(sum([x[idx] for x in posesL]))
		for idx in range(distR):
			dp_listR.append(sum([x[idx] for x in posesR]))
		tdpL=median(dp_listL)
		tdpR=median(dp_listR)
		return tdpL, tdpR
	
	tdp1L, tdp1R = subr(chr1, pos1, tbam_file)
	ndp1L, ndp1R = subr(chr1, pos1, nbam_file)
	
	if ter1=='5':
		depth_ratio_bp1_in = (tdp1R+1)/(ndp1R+1)
		depth_ratio_bp1_out = (tdp1L+1)/(ndp1L+1)
	elif ter1=='3':
		depth_ratio_bp1_in = (tdp1L+1)/(ndp1L+1)
		depth_ratio_bp1_out = (tdp1R+1)/(ndp1R+1)
		
	tdp2L, tdp2R = subr(chr2, pos2, tbam_file)
	ndp2L, ndp2R = subr(chr2, pos2, nbam_file)
	
	if ter2=='5':
		depth_ratio_bp2_in = (tdp2R+1)/(ndp2R+1)
		depth_ratio_bp2_out = (tdp2L+1)/(ndp2L+1)
	elif ter2=='3':
		depth_ratio_bp2_in = (tdp2L+1)/(ndp2L+1)
		depth_ratio_bp2_out = (tdp2R+1)/(ndp2R+1)
		
	depth_ratio_change_bp1 = depth_ratio_bp1_in - depth_ratio_bp1_out
	depth_ratio_change_bp2 = depth_ratio_bp2_in - depth_ratio_bp2_out
	
	return depth_ratio_change_bp1, depth_ratio_change_bp2
	
	
def collect_fetchrange(chr1, pos1, ter1, chr2, pos2, ter2, sv_type, dist, chr_size):
	def merge_ranges(rangelist):
		chrom = rangelist[0][0]
		start = max( 0, min([x[1] for x in rangelist]) )
		end = min( max([x[2] for x in rangelist]), chr_size[chrom] )
		return chrom, start, end
		
		
	fetchrange_list_tumor_chr1 = list()
	fetchrange_list_tumor_chr2 = list()
	fetchrange_list_normal_chr1 = list()
	fetchrange_list_normal_chr2 = list()
	
	# step1
	wrapper_step1 = step1_wrapper(chr1=chr1, pos1=pos1, ter1=ter1, chr2=chr2, pos2=pos2, ter2=ter2, 
									 sv_type=sv_type, dist=dist)
	wrapper_step1.get_fetchrange()
	
	fetchrange_list_tumor_chr1.append(wrapper_step1.fetchrange_dict['chr1']) # chr1, start1 - 1, end1
	fetchrange_list_normal_chr1.append(wrapper_step1.fetchrange_dict['chr1'])
	fetchrange_list_tumor_chr2.append(wrapper_step1.fetchrange_dict['chr2']) # chr2, start2 - 1, end2
	fetchrange_list_normal_chr2.append(wrapper_step1.fetchrange_dict['chr2'])
	
	# step3
	# calc_final_count
	wrapper_calcFinalCount = calc_final_count(chr1=chr1, pos1=pos1, ter1=ter1, chr2=chr2, pos2=pos2, ter2=ter2)
	wrapper_calcFinalCount.get_fetchrange()
	
	fetchrange_list_tumor_chr1.append(wrapper_calcFinalCount.fetchrange_dict['chr1']) # chr1, pos1_start - 1, pos1_end
	fetchrange_list_normal_chr1.append(wrapper_calcFinalCount.fetchrange_dict['chr1'])
	fetchrange_list_tumor_chr2.append(wrapper_calcFinalCount.fetchrange_dict['chr2']) # chr2, pos2_start - 1, pos2_end
	fetchrange_list_normal_chr2.append(wrapper_calcFinalCount.fetchrange_dict['chr2'])
	
	# count_frag_num
	wrapper_countFragNum = count_frag_num(chr1=chr1, pos1=pos1)
	wrapper_countFragNum.get_fetchrange()
	fetchrange_list_normal_chr1.append(wrapper_countFragNum.fetchrange) # chr1, pos1 - 1, pos1
	wrapper_countFragNum = count_frag_num(chr1=chr2, pos1=pos2)
	wrapper_countFragNum.get_fetchrange()
	fetchrange_list_normal_chr2.append(wrapper_countFragNum.fetchrange) # chr2, pos2 - 1, pos2
	
	# step4
	wrapper_findPnsc = find_pnsc(chr1=chr1, pos1=pos1, ter1=ter1, chr2=chr2, pos2=pos2, ter2=ter2)
	wrapper_findPnsc.get_fetchrange()
	fetchrange_list_tumor_chr1.append(wrapper_findPnsc.fetchrange_dict['tumor']) # chr1, pos1_start - 1, pos1_end
	fetchrange_list_normal_chr1.append(wrapper_findPnsc.fetchrange_dict['normal']) # chr1, max(1, pos1-1-n_sr), pos1+n_sr
	
	wrapper_findPnsc = find_pnsc(chr1=chr2, pos1=pos2, ter1=ter2, chr2=chr1, pos2=pos1, ter2=ter1)
	wrapper_findPnsc.get_fetchrange()
	fetchrange_list_tumor_chr2.append(wrapper_findPnsc.fetchrange_dict['tumor']) # chr2, ...
	fetchrange_list_normal_chr2.append(wrapper_findPnsc.fetchrange_dict['normal']) # chr2, 
	
	# step5
	wrapper_amountDiscordant = amount_discordant(chr1=chr1, pos1=pos1)
	wrapper_amountDiscordant.get_fetchrange()
	fetchrange_list_tumor_chr1.append(wrapper_amountDiscordant.fetchrange)
	fetchrange_list_normal_chr1.append(wrapper_amountDiscordant.fetchrange)
	
	wrapper_amountDiscordant = amount_discordant(chr1=chr2, pos1=pos2)
	wrapper_amountDiscordant.get_fetchrange()
	fetchrange_list_tumor_chr2.append(wrapper_amountDiscordant.fetchrange)
	fetchrange_list_normal_chr2.append(wrapper_amountDiscordant.fetchrange)
	
	# merge partial ranges
	result = dict()
	result['tumor'] = dict() ; result['normal'] = dict()
	
	result['tumor']['chr1'] = merge_ranges(fetchrange_list_tumor_chr1)
	result['tumor']['chr2'] = merge_ranges(fetchrange_list_tumor_chr2)
	result['normal']['chr1'] = merge_ranges(fetchrange_list_normal_chr1)
	result['normal']['chr2'] = merge_ranges(fetchrange_list_normal_chr2)
	
	return result


def get_rplist(fetchrange_collection, tbam_file, nbam_file, ter1, r_limit_dp, chr_list, primary_contigs, nonprim, readlength = 150):
	
	def filter_fun(read):
		return not (
			read.cigarstring == None or
			read.is_unmapped or
			read.is_duplicate
		)
	
	def selective_fetch(bam, fetchrange, cutoff):
		if cutoff == None:
			result = list()
			readcount = None
			for read in bam.fetch(*fetchrange):
				if filter_fun(read):
					result.append(read)
		else:
			result = list()
			readcount = bam.count(*fetchrange, read_callback = filter_fun)
			if readcount > cutoff:
				prop = cutoff / readcount
				for read in bam.fetch(*fetchrange):
					if filter_fun(read):
						if random.random() < prop:
							result.append(read)
			else:
				for read in bam.fetch(*fetchrange):
					if filter_fun(read):
						result.append(read)
					
		return result, readcount
	
	bamdict = dict()
	bamdict['tumor'] = tbam_file ; bamdict['normal'] = nbam_file
	
	# fetch reads with random choice
	readlist = dict()
	readlist['tumor'] = dict() ; readlist['normal'] = dict()
	for bamtype in readlist:
		for chromtype in ('chr1', 'chr2'):
			fetchrange_width = fetchrange_collection[bamtype][chromtype][2] - fetchrange_collection[bamtype][chromtype][1]
			if r_limit_dp == None:
				cutoff = None
			else:
				cutoff = int( r_limit_dp * math.ceil(fetchrange_width/readlength) )
				
			readlist[bamtype][chromtype], readcount = \
				selective_fetch(
					bamdict[bamtype], 
					fetchrange_collection[bamtype][chromtype],
					cutoff
				)
				
	# make readplus objects
	rplist_dict = dict()
	rplist_dict['tumor'] = dict() ; rplist_dict['normal'] = dict()
	for bamtype in rplist_dict:
		for chromtype in ('chr1', 'chr2'):
			rplist_dict[bamtype][chromtype] = [ readplus(read, ter1, chr_list, primary_contigs, nonprim) for read in readlist[bamtype][chromtype] ]
	
	return rplist_dict
	

def variant_filter(chr1, pos1, chr2, pos2, pr_satellite, primary_contigs, sv_type, nonprim, chr_size): # Filter out if True is returned
	# skip if overlapping with satellite
	pr_variant_chr1 = pr.PyRanges(chromosomes = chr1, starts = [pos1 - 1], ends = [pos1])
	pr_variant_chr2 = pr.PyRanges(chromosomes = chr2, starts = [pos2 - 1], ends = [pos2])
	if (not pr_satellite.intersect(pr_variant_chr1).empty) or (not pr_satellite.intersect(pr_variant_chr2).empty):
		printErr('satellite')
		return True

	# skip if any of the two breakends are not in the primary contigs
	if not nonprim:
		if (chr1 not in primary_contigs) or (chr2 not in primary_contigs):
			printErr('non-primary contig')
			return True

	# skip if sv_type is INS
	if sv_type == 'INS': 
		return True

	# skip if pos is greater than the length of chr
	if pos1 > chr_size[chr1] or pos2 > chr_size[chr2]:
		printErr('"pos" out of chromosome length range')
		return True


	return False
