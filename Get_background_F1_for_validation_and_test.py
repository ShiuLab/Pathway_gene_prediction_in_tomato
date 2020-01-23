import sys,os
import numpy as np
def Fmeasure(TP,FP,FN):
	if TP+FP != 0:
		Pre = float(TP)/(TP+FP)
	if TP+FN != 0:
		Rec = float(TP)/(TP+FN)
	if TP+FP != 0 and TP+FN != 0 and Pre != 0:
		F1 = (2*Pre*Rec)/float(Pre+Rec)
	else: 
		F1 = 0
	return(F1)

F = {}
for i in range(1,6):
	inp = open('Pathway_genes_for_5_validation_set%s.txt'%i,'r').readlines()
	R = {}
	for inl in inp:
		pathway = inl.split('\t')[0]
		if pathway not in R:
			R[pathway] = []
		R[pathway].append(inl.split('\t')[1].strip())
		if pathway not in F:
			F[pathway] = []
	for pathway in R:
		TP = len(R[pathway]) * float(len(R[pathway]))/(len(inp)-len(R[pathway]))
		FN = len(R[pathway]) * (1-float(len(R[pathway]))/(len(inp)-len(R[pathway])))
		FP = (len(inp)-len(R[pathway])) * float(len(R[pathway]))/(len(inp)-len(R[pathway]))
		F[pathway].append(Fmeasure(TP,FP,FN))
		
out = open('Background_F1_for_validation.txt','w')
for pathway in F:
	out.write('%s\t%s\n'%(pathway,np.mean(F[pathway])))
	
out.close()
		
inp = open('Pathway_genes_for_testing.txt','r').readlines()
R = {}
F = {}
for inl in inp:
	pathway = inl.split('\t')[0]
	if pathway not in R:
		R[pathway] = []
	R[pathway].append(inl.split('\t')[1].strip())
	if pathway not in F:
		F[pathway] = []
for pathway in R:
	TP = len(R[pathway]) * float(len(R[pathway]))/(len(inp)-len(R[pathway]))
	FN = len(R[pathway]) * (1-float(len(R[pathway]))/(len(inp)-len(R[pathway])))
	FP = (len(inp)-len(R[pathway])) * float(len(R[pathway]))/(len(inp)-len(R[pathway]))
	F[pathway].append(Fmeasure(TP,FP,FN))

out = open('Background_F1_for_testing.txt','w')
for pathway in F:
	out.write('%s\t%s\n'%(pathway,F[pathway][0]))
	
out.close()

inp = open('Sly_pathway_annotation_20190117_with_expression_5_members_nonoverlapping.txt','r').readlines()
R = {}
F = {}
for inl in inp:
	pathway = inl.split('\t')[0]
	if pathway not in R:
		R[pathway] = []
	R[pathway].append(inl.split('\t')[1].strip())
	if pathway not in F:
		F[pathway] = []
for pathway in R:
	TP = len(R[pathway]) * float(len(R[pathway]))/(len(inp)-len(R[pathway]))
	FN = len(R[pathway]) * (1-float(len(R[pathway]))/(len(inp)-len(R[pathway])))
	FP = (len(inp)-len(R[pathway])) * float(len(R[pathway]))/(len(inp)-len(R[pathway]))
	F[pathway].append(Fmeasure(TP,FP,FN))

out = open('Background_F1_for_all_nonoverlapping.txt','w')
for pathway in F:
	out.write('%s\t%s\n'%(pathway,F[pathway][0]))
	
out.close()

