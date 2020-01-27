import sys,os
import pandas as pd
import numpy as np
def Fmeasure(TP,FP,FN):
	if TP+FP != 0:
		Pre = float(TP)/(TP+FP)
	if TP+FN != 0:
		Rec = float(TP)/(TP+FN)
	if TP+FP != 0 and TP+FN != 0 and Pre != 0 and TP != 0:
		F1 = (2*Pre*Rec)/float(Pre+Rec)
	else: 
		F1 = 0
	return(F1)

path = sys.argv[1]
expression = sys.argv[2]
Result = {}
Gene_prediction = {}  # Gene_prediction[gene1][Median_1] = pathway
for cv_number in range(1,6):
	gene_list = open('Pathway_genes_for_5_validation_set%s.txt'%cv_number,'r').readlines()
	D = {} #D[gene] = pathway
	P = {} #P[pathway] = [gene1,gene2,...]
	for inl in gene_list:
		gene = inl.split('\t')[1].strip()
		pathway = inl.split('\t')[0]
		if gene not in D:
			D[gene] = pathway
		else:
			print('Redundance of %s'%gene)
		if pathway not in P:
			P[pathway] = []
		P[pathway].append(gene)
	training = open('Pathway_genes_for_5_training_set%s.txt'%cv_number,'r').readlines()
	D_training = {} #D_training[gene] = pathway
	P_training = {} #P_training[pathway] = [gene1,gene2,...]
	for inl in training:
		gene = inl.split('\t')[1].strip()
		pathway = inl.split('\t')[0]
		if gene not in D_training:
			D_training[gene] = pathway
		else:
			print('Redundance of %s'%gene)
		if pathway not in P_training:
			P_training[pathway] = []
		P_training[pathway].append(gene)
	df = pd.read_csv(path + expression, sep='\t', index_col = 0, header = 0)
	if expression.startswith('Multiclass_MR'):
		df = df.iloc[:,1:2172]
		colname = df.columns.tolist()
		colname = [c.replace(expression.split('Multiclass_')[1].split('complete.txt')[0],'') for c in colname]
		df.columns = colname
	rownames = df.index.tolist()
	Predicted_max = {} ### Predicted_max[pathwayA] = [gene1,gene2]
	Predicted_median = {} ### Predicted_median[pathwayA] = [gene1,gene2]
	for gene1 in D:
		Median = {}
		Max = {}
		for pathway in P_training:
			PCC = []
			for gene2 in P_training[pathway]:
				if gene1 != gene2:
					try:
						PCC.append(float(df.loc[gene1,gene2]))  ### PCC values of gene1 to all genes in a pathway
						if float(df.loc[gene1,gene2]) not in Max:
							Max[float(df.loc[gene1,gene2])] = []  
						Max[float(df.loc[gene1,gene2])].append(gene2) ### unique PCC values of gene1 to all other genes
					except:
						print('Na for correlation between %s and %s'%(gene1,gene2))
			try:
				if np.median(PCC) not in Median:
					Median[np.median(PCC)] = []
				Median[np.median(PCC)].append(pathway) ### unique Median PCC values of gene1 to all pathways
			except:
				print('Na for correlation between %s'%(gene1))
		try:
			perdicted_as_max = Max[max(Max.keys())]
			if len(perdicted_as_max)==1:
				pathway_predicted = D_training[perdicted_as_max[0]]
				if pathway_predicted not in Predicted_max:
					Predicted_max[pathway_predicted] = []
				if gene1 not in Predicted_max[pathway_predicted]:
					Predicted_max[pathway_predicted].append(gene1)
				if gene1 not in Gene_prediction:
					Gene_prediction[gene1] = {}
				if 'Max_%s'%cv_number not in Gene_prediction[gene1]:
					Gene_prediction[gene1]['Max_%s'%cv_number] = pathway_predicted
		except:
			print('No assignment for %s'%gene1)
		try:
			perdicted_as_median = Median[max(Median.keys())]
			if len(perdicted_as_median) == 1:
				if perdicted_as_median[0] not in Predicted_median:
					Predicted_median[perdicted_as_median[0]] = []
				if gene1 not in Predicted_median[perdicted_as_median[0]]:
					Predicted_median[perdicted_as_median[0]].append(gene1)
				if gene1 not in Gene_prediction:
					Gene_prediction[gene1] = {}
				if 'Median_%s'%cv_number not in Gene_prediction[gene1]:
					Gene_prediction[gene1]['Median_%s'%cv_number] = perdicted_as_median[0]
		except:
			print('No assignment for %s'%gene1)
	for pathway in P:
		if pathway not in Result:
			Result[pathway] = {}
		if 'Median' not in Result[pathway]:
			Result[pathway]['Median'] = []
		if 'Max' not in Result[pathway]:
			Result[pathway]['Max'] = []
		TP_median = 0
		FP_median = 0
		FN_median = 0
		TN_median = 0
		TP_max = 0
		FP_max = 0
		FN_max = 0
		TN_max = 0
		if pathway in Predicted_median:
			for gene in Predicted_median[pathway]:
				if gene in P[pathway]:
					TP_median += 1
				if gene not in P[pathway]:
					FP_median += 1
			FN_median = len(P[pathway]) - TP_median
			TN_median = len(D) - TP_median - FP_median - FN_median
		else:
			TP_median = 0
			FP_median = 0
			FN_median = len(P[pathway]) - TP_median
			TN_median = len(D) - TP_median - FP_median - FN_median
		if pathway in Predicted_max:
			for gene in Predicted_max[pathway]:
				if gene in P[pathway]:
					TP_max += 1
				if gene not in P[pathway]:
					FP_max += 1
			FN_max = len(P[pathway]) - TP_max
			TN_max = len(D) - TP_max - FP_max - FN_max
		else:
			TP_max = 0
			FP_max = 0
			FN_max = len(P[pathway]) - TP_max
			TN_max = len(D) - TP_max - FP_max - FN_max
		F1_max = Fmeasure(TP_max,FP_max,FN_max)
		F1_median = Fmeasure(TP_median,FP_median,FN_median)
		Result[pathway]['Max'].append(F1_max)
		Result[pathway]['Median'].append(F1_median)

out = open('./Result_with_geen_prediction/'+ expression.split('.txt')[0].split('/')[-1] + '_naive_median_max_validation.txt','w')
out.write('Pathway\tMedian_F1_average\tMedian_F1_std\tMax_F1_average\tMax_F1_std\n')
for pathway in Result:
	out.write('%s\t%s\t%s\t%s\t%s\n'%(pathway,np.mean(Result[pathway]['Median']),np.std(Result[pathway]['Median']),np.mean(Result[pathway]['Max']),np.std(Result[pathway]['Max'])))
			
out.close()	

inp = open('Sly_pathway_annotation_20190117_with_expression_5_members_nonoverlapping.txt','r').readlines()
for inl in inp:
	gene = inl.split('\t')[1].strip()
	pathway = inl.split('\t')[0]
	if gene in Gene_prediction:
		Gene_prediction[gene]['Annotated'] = pathway

	
gene_prediction = pd.DataFrame.from_dict(Gene_prediction, orient='index')  ### make a dataframe from a dictionary
gene_prediction.to_csv('./Result_with_geen_prediction/'+ expression.split('.txt')[0] + '_gene_predicted.txt', index=True, header=True,sep="\t")
	
			
