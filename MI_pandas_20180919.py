'''
export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH
### should ssh dev-intel16-k80

input1: expression data, Fold change or FPKM
input2: start
input3: stop
'''
import sys,os,argparse
import pandas as pd
import numpy as np
import math
from sklearn.metrics.cluster import normalized_mutual_info_score

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn


def main():
	parser = argparse.ArgumentParser(description='This code is for calculating the Mutual information for a gene expression matrix. Because it is too slow to calculate a big matrix, you may want to do the job for a subset of your data, defined by the starting row number and the stoping row number.')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-file', help='Expression matrix', required=True)
	req_group.add_argument('-path', help='path to the Expression matrix', required=True)
	req_group.add_argument('-start', help='where the subset starts', required=True)
	req_group.add_argument('-stop', help='where the subset stops', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	

	file = args.file
	start = int(args.start)
	stop = int(args.stop)

	out = open('MI_%s_%s_%s'%(file,start,stop),'w')

	df = pd.read_csv(path+file, sep='\t', index_col = 0, header = 0)
	D = {} ###
	rowname = df.index.tolist()
	title = 'gene'
	for name in rowname:
		title = title + '\t' + name
	out.write(title + '\n')
	out.flush()

	x = start -1
	while x < stop:
		gene1 = rowname[x]
		result = gene1
		for gene2 in rowname:
			MI = float(normalized_mutual_info_score(df.loc[gene1,:],df.loc[gene2,:]))
			result = result + '\t%s'%MI
		out.write(result + '\n')
		out.flush()
		x += 1

	out.close()

if __name__ == '__main__':
	main()

