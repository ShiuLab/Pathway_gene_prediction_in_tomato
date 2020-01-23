import sys,os
import pandas as pd
type = sys.argv[1]
D = {}
for file in os.listdir('./'):
	if file.endswith('results.txt'):
		inp = open(file,'r').readlines()
		D[file.split('_RF_results.txt')[0]] = {}
		x = 0
		while x < len(inp) and not inp[x].startswith('HO F1_macro'):
			x += 1
				
		x += 1
		while x < len(inp):
			inl = inp[x]
			if inl.split('\t')[0].endswith('_F1'):
				D[file.split('_RF_results.txt')[0]][inl.split('\t')[0].split('_F1')[0]] = float(inl.split('\t')[1])
			x += 1
				
res = pd.DataFrame.from_dict(D)
res.to_csv('../F1_matrix_on_test_pathways_%s.txt'%type, index=True, header=True,sep="\t")