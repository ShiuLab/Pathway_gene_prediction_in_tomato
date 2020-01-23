import sys,os
import random
import copy
from random import shuffle
inp = open('Sly_pathway_annotation_20190117_with_expression_5_members_nonoverlapping.txt','r').readlines()
P = {}  ###P[pathway] = [gene1,gene2,...]
for inl in inp:
	pa = inl.strip().split('\t')[0]
	gene = inl.split('\t')[1].strip()
	if pa not in P:
		P[pa] = []
	P[pa].append(gene)

size = open('Nonoverlapping_pathway_size.txt','w')
for pa in P:
	size.write('%s\t%s\n'%(pa,len(P[pa])))
	
size.close()
	
test = open('Pathway_genes_for_testing.txt','w')
T = {}
P_copy = copy.deepcopy(P)
for pa in P:
	if len(P[pa])>25:
		random_t = random.sample(P[pa],5)
		for t in random_t:
			if pa not in T:
				T[pa] = []
			T[pa].append(t)
			P_copy[pa].remove(t)
			test.write('%s\t%s\n'%(pa,t))

test.close()

training1 = open('Pathway_genes_for_5_training_set1.txt','w')
validation1 = open('Pathway_genes_for_5_validation_set1.txt','w')
training2 = open('Pathway_genes_for_5_training_set2.txt','w')
validation2 = open('Pathway_genes_for_5_validation_set2.txt','w')
training3 = open('Pathway_genes_for_5_training_set3.txt','w')
validation3 = open('Pathway_genes_for_5_validation_set3.txt','w')
training4 = open('Pathway_genes_for_5_training_set4.txt','w')
validation4 = open('Pathway_genes_for_5_validation_set4.txt','w')
training5 = open('Pathway_genes_for_5_training_set5.txt','w')
validation5 = open('Pathway_genes_for_5_validation_set5.txt','w')
for pa in P_copy:
	gene_list = list(P_copy[pa])
	shuffle(gene_list)
	residue = len(gene_list)%5
	aa = [0,0,0,0,0]
	for i in range(5-residue,5):
		aa[i] += 1
	for n in range(0,int(len(gene_list)/5)+aa[0]):
		validation1.write('%s\t%s\n'%(pa,gene_list[n]))
		training2.write('%s\t%s\n'%(pa,gene_list[n]))
		training3.write('%s\t%s\n'%(pa,gene_list[n]))
		training4.write('%s\t%s\n'%(pa,gene_list[n]))
		training5.write('%s\t%s\n'%(pa,gene_list[n]))
	for n in range(int(len(gene_list)/5)+aa[0],int(len(gene_list)/5) * 2 +aa[0]+aa[1]):
		validation2.write('%s\t%s\n'%(pa,gene_list[n]))
		training1.write('%s\t%s\n'%(pa,gene_list[n]))
		training3.write('%s\t%s\n'%(pa,gene_list[n]))
		training4.write('%s\t%s\n'%(pa,gene_list[n]))
		training5.write('%s\t%s\n'%(pa,gene_list[n]))
	for n in range(int(len(gene_list)/5) * 2 +aa[0]+aa[1] ,int(len(gene_list)/5) * 3 +aa[0]+aa[1]+aa[2]):
		validation3.write('%s\t%s\n'%(pa,gene_list[n]))
		training1.write('%s\t%s\n'%(pa,gene_list[n]))
		training2.write('%s\t%s\n'%(pa,gene_list[n]))
		training4.write('%s\t%s\n'%(pa,gene_list[n]))
		training5.write('%s\t%s\n'%(pa,gene_list[n]))
	for n in range(int(len(gene_list)/5) * 3 +aa[0]+aa[1]+aa[2] ,int(len(gene_list)/5) * 4 +aa[0]+aa[1]+aa[2]+aa[3]):
		validation4.write('%s\t%s\n'%(pa,gene_list[n]))
		training1.write('%s\t%s\n'%(pa,gene_list[n]))
		training2.write('%s\t%s\n'%(pa,gene_list[n]))
		training3.write('%s\t%s\n'%(pa,gene_list[n]))
		training5.write('%s\t%s\n'%(pa,gene_list[n]))
	for n in range(int(len(gene_list)/5) * 4 +aa[0]+aa[1]+aa[2]+aa[3],int(len(gene_list)/5) * 5 +aa[0]+aa[1]+aa[2]+aa[3]+aa[4]):
		validation5.write('%s\t%s\n'%(pa,gene_list[n]))
		training1.write('%s\t%s\n'%(pa,gene_list[n]))
		training2.write('%s\t%s\n'%(pa,gene_list[n]))
		training3.write('%s\t%s\n'%(pa,gene_list[n]))
		training4.write('%s\t%s\n'%(pa,gene_list[n]))	
	
validation1.close()
training1.close()
validation2.close()
training2.close()
validation3.close()
training3.close()
validation4.close()
training4.close()
validation5.close()
training5.close()







	