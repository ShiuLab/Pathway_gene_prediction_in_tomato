'''
export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH
### should ssh dev-intel16-k80
'''
import sys,os,argparse
import pandas as pd
import numpy as np

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for calculating the Spearman ρ for a gene expression matrix')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-file', help='Expression matrix', required=True)
	req_group.add_argument('-path', help='path to the Expression matrix', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	

	file = args.file
	SAVE = 'Spearman_' + file

	df = pd.read_csv(path+file, sep='\t', index_col = 0, header = 0)  
	pcc = df.T.corr(method='spearman').round(3)
	pcc.to_csv(SAVE, index=True, header=True,sep="\t")

if __name__ == '__main__':
	main()

