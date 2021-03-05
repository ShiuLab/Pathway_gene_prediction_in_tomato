import sys,os,argparse
import numpy as np

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for calculating the Mutual Ranks for gene pairs')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-pathway_anotation', help='pathway annotation', required=True)
	req_group.add_argument('-exp', help='Expression matrix', required=True)
	req_group.add_argument('-method', help='methods for coexpression. pcc, spearman, MI, or corpcor', required=True)

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

	Tem = open('Multiclass_%s_%s.txt'%(method,S[exp]),'r').readlines()
	out = open('Multiclass_MR_%s_%s.txt'%(method,S[exp]),'w')
	title = 'Gene\tPathway'
	title_tem = Tem[0].strip().split('\t')
	for t in title_tem[2:]:
		title = title + '\tMR_%s'%t
	out.write(title + '\n')
	D_R = {}  ###D_R[gene][target] = [Rank]
	for inl in Tem[1:]:  
		gene = inl.split('\t')[0]
		pathway = inl.split('\t')[1]
		res_tem = '\t'.join(inl.split('\t')[0:len(G)+2])
		D = {}  	  
		D_R[gene] = {}
		D_1 = []
		result_tem = res_tem.split('\t')
		for i in range(2,len(G)+2):	
			g = title_tem[i].split('_')[-2] + '_' + title_tem[i].split('_')[-1]
			D[g] = float(result_tem[i])
			if gene != g:
				D_1.append(float(result_tem[i]))
		D_2 = np.array(D_1)
		for i in range(2,len(G)+2):		
			g = title_tem[i].split('_')[-2] + '_' + title_tem[i].split('_')[-1]
			if gene != g:
				rank1_2 = sum((D_2>D[g]).astype(np.int)) + 1
				D_R[gene][g] = rank1_2
			else:
				D_R[gene][g] = 0
			
	for inl in Tem[1:]:  
		gene = inl.split('\t')[0]
		pathway = inl.split('\t')[1]
		result = '%s\t%s'%(gene,pathway)
		R = {}
		for i in range(2,len(G)+2):	
			g = title_tem[i].split('_')[-2] + '_' + title_tem[i].split('_')[-1]
			value = (D_R[gene][g]*D_R[g][gene])**(0.5)
			result = result + '\t%s'%value	
		for i in range(len(G)+2,len(title_tem),2):	
			p = title_tem[i].split('_')[-2]
			R[p] = []
			for g in P[p]:
				if g != gene:
					value = (D_R[gene][g]*D_R[g][gene])**(0.5)
					R[p].append(float(value))				
		for i in range(len(G)+2,len(title_tem),2):	
			p = title_tem[i].split('_')[-2]
			max = np.max(R[p])
			median = np.median(R[p])
			result = result + '\t%s\t%s'%(max,median)
		out.write(result + '\n')
		out.flush()

	out.close()


if __name__ == '__main__':
	main()
