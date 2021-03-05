import sys,os,argparse
import os,sys
import pandas as pd
import numpy as np
import random
from scipy import stats
import itertools
import math
import sklearn
from sklearn.metrics.cluster import normalized_mutual_info_score
def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for calculating the gene-to-gene and gene-to-pathway expression similarities. if the method is pcc, spearman, MI, the original expression matrix will be loaded, e.g., FPKM or TPM matrix. If the method is corpcor, then the corpcor gene-to-gene matrix will be loaded, thus make sure you run the script corpcor.r first to get the partial coexpression matrix.')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-pathway_anotation', help='pathway annotation', required=True)
	req_group.add_argument('-exp', help='Expression matrix', required=True)
	req_group.add_argument('-method', help='methods for coexpression. pcc, spearman, MI or corpcor', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	

	pa = args.pathway_anotation
	exp = args.exp
	method = args.method
	sample = open('Short_name_for_expression_data.txt','r').readlines()
	S = {}  ### short name for expression data file
	for inl in sample:
		S[inl.split('\t')[0]] = inl.strip().split('\t')[1]

	if method in ['pcc','spearman','MI']:
		df = pd.read_csv(exp, sep='\t', index_col = 0, header = 0)
	if method == 'corpcor':
		df = pd.read_csv('corpcor_' + exp, sep='\t', index_col = 0, header = 0)

	inp = open(pa,'r').readlines()
	P = {} ### pathway
	G = {}
	for inl in inp:
		pathway = inl.split('\t')[0]
		gene = inl.split('\t')[1].strip()
		if pathway not in P:
			P[pathway] = []
		P[pathway].append(gene)
		if gene not in G:
			G[gene] = []
		G[gene].append(pathway)

	out = open('Multiclass_%s_%s.txt'%(method,S[exp]),'w')
	title = 'Gene\tPathway'
	for gene in sorted(G.keys()):
		title = '%s\t%s_%s_%s'%(title,method,S[exp],gene)
	for pathway in sorted(P.keys()):
		title = '%s\t%s_%s_%s_max\t%s_%s_%s_median'%(title,method,S[exp],pathway,method,S[exp],pathway)
	out.write(title + '\n')
	title_tem = title.split('\t')

	for gene in sorted(G.keys()):
		result = '%s\t%s'%(gene,G[gene][0])
		for i in range(2,(len(G)+2)):
			target = title_tem[i].split('_')[-2] + '_' + title_tem[i].split('_')[-1]
			if method == 'pcc':
				value = stats.pearsonr(df.loc[gene,:],df.loc[target,:])[0]
			if method == 'spearman':
				value = stats.spearmanr(df.loc[gene,:],df.loc[target,:])[0]
			if method == 'MI':
				value = normalized_mutual_info_score(df.loc[gene,:],df.loc[target,:])
			if method == 'corpcor':
				value = df.loc[gene,target]
			if math.isnan(value):
				value = 0
			result = result + '\t%s'%value
		D = {}
		result_tem = result.split('\t')
		for i in range(2,(len(G)+2)):		
			D[title_tem[i].split('_')[-2] + '_' + title_tem[i].split('_')[-1]] = result_tem[i]
		R = {}
		for pathway in sorted(P.keys()):
			R[pathway] = []
			for g in P[pathway]:
				if g != gene:
					R[pathway].append(float(D[g]))
		for pathway in sorted(P.keys()):
			max = np.max(R[pathway])
			median = np.median(R[pathway])
			result = result + '\t%s\t%s'%(max,median)
		out.write(result + '\n')
		out.flush()
			
	out.close()

if __name__ == '__main__':
	main()
