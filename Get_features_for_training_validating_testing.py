import sys,os
import pandas as pd
import numpy as np
path = sys.argv[1] # the pathway where the co-expression matrix are
files = sys.argv[2] # the co-expression matrix
number = sys.argv[3] # the numbers for cross-validation
expression = pd.read_csv(path + files, sep='\t', index_col = 0, header = 0)
expression = expression.drop(['Pathway'], axis=1)
shortname = '_'.join(expression.columns[0].split('_')[0:-2]) + '_'
expression.columns = expression.columns.str.replace(shortname,"")

training = open('Pathway_genes_for_5_training_set%s.txt'%number,'r').readlines()
validation = open('Pathway_genes_for_5_validation_set%s.txt'%number,'r').readlines()
test = open('Pathway_genes_for_testing.txt','r').readlines()
P_train = {}
P_validation = {}
P_test = {}
for inl in training:
	tem = inl.strip().split('\t')
	if tem[0] not in P_train:
		P_train[tem[0]] = []
	P_train[tem[0]].append(tem[1])
	
for inl in validation:
	tem = inl.strip().split('\t')
	if tem[0] not in P_validation:
		P_validation[tem[0]] = []
	P_validation[tem[0]].append(tem[1])
	
for inl in test:
	tem = inl.strip().split('\t')
	if tem[0] not in P_test:
		P_test[tem[0]] = []
	P_test[tem[0]].append(tem[1])
	
pathway = sorted(P_train.keys())
out = open(shortname+'CV_%s_features.txt'%number,'w')
title = 'Gene\tClass'
for p in pathway:
	title = title + '\t%s%s_median\t%s%s_max'%(shortname,p, shortname,p)
	
out.write(title + '\n')

for p in pathway:
	for gene1 in P_train[p]:
		res = '%s\t%s'%(gene1,p)
		for p2 in pathway:
			value = []
			for gene2 in P_train[p2]:
				if gene1!=gene2:
					value.append(round(expression.loc[gene1,gene2],3))
			if 'MR' in files:
				res = res + '\t%s\t%s'%(round(np.median(value),3),round(np.min(value),3))
			else:
				res = res + '\t%s\t%s'%(round(np.median(value),3),round(np.max(value),3))
		out.write(res + '\n')
				
for p in pathway:
	for gene1 in P_validation[p]:
		res = '%s\t%s'%(gene1,p)
		for p2 in pathway:
			value = []
			for gene2 in P_train[p2]:
				value.append(round(expression.loc[gene1,gene2],3))
			if 'MR' in files:
				res = res + '\t%s\t%s'%(round(np.median(value),3),round(np.min(value),3))
			else:
				res = res + '\t%s\t%s'%(round(np.median(value),3),round(np.max(value),3))
		out.write(res + '\n')
					
for p in P_test.keys():
	for gene1 in P_test[p]:
		res = '%s\t%s'%(gene1,p)
		for p2 in pathway:
			value = []
			for gene2 in P_train[p2]:
				value.append(round(expression.loc[gene1,gene2],3))
			if 'MR' in files:
				res = res + '\t%s\t%s'%(round(np.median(value),3),round(np.min(value),3))
			else:
				res = res + '\t%s\t%s'%(round(np.median(value),3),round(np.max(value),3))
		out.write(res + '\n')

out.close()		
