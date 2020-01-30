import sys, os
import pandas as pd
import numpy as np
from datetime import datetime
import time
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import cross_val_predict
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

def Performance_MC(y, pred):
	from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
	pred = pd.DataFrame(pred,index=y.index)
	df_tem = pd.concat([y,pred],axis=1)
	df_tem.columns = ['y_true','y_pred']
	f1 = []
	for y_list in y.unique():
		P = df_tem[df_tem.y_true==y_list]
		TP = P[P.y_pred == y_list]
		FN = P.shape[0] - TP.shape[0]
		N = df_tem[df_tem.y_true!=y_list]
		FP = N[N.y_pred == y_list]
		f1_tem = Fmeasure(TP.shape[0],FP.shape[0],FN)
		f1.append(f1_tem)
	d = {'F1': f1}
	result = pd.DataFrame(data=d)
	result.index = y.unique()
	return(result)

n = 0
for files in os.listdir('./'):
	if files.endswith('_complete_gene_predicted_test.txt'):
		name = files.split('_complete_gene_predicted_test.txt')[0].split('Multiclass_')[1]
		df = pd.read_csv(files,index_col=0,header=0,sep='\t')
		df2 = pd.DataFrame([['PWY', 'PWY', 'PWY', 'PWY', 'PWY', 'PWY', 'PWY', 'PWY', 'PWY', 'PWY', 'PWY']],columns=list(df.columns.tolist()))
		df = df.append(df2,ignore_index=True)
		df = df.fillna('PWY')
		y = df.Annotated
		m_max = 0
		for pred in [df.Max_1,df.Max_2,df.Max_3,df.Max_4,df.Max_5]:
			a = Performance_MC(y,pred)
			if m_max == 0:
				max_result = a
			else:
				max_result = pd.concat([max_result,a],ignore_index=True,axis=1)
			m_max += 1
		max_result['Mean_F1'] = max_result.mean(axis=1)
		d = {name:max_result['Mean_F1']}
		dd = pd.DataFrame(data=d)
		if n == 0:
			F1_max = dd
		else:
			F1_max = pd.concat([F1_max,dd],axis=1)
		m_median = 0
		for pred in [df.Median_1,df.Median_2,df.Median_3,df.Median_4,df.Median_5]:
			a = Performance_MC(y,pred)
			if m_median == 0:
				median_result = a
			else:
				median_result = pd.concat([median_result,a],ignore_index=True,axis=1)
			m_median += 1
		median_result['Mean_F1'] = median_result.mean(axis=1)
		d = {name:median_result['Mean_F1']}
		dd = pd.DataFrame(data=d)
		if n == 0:
			F1_median = dd
		else:
			F1_median = pd.concat([F1_median,dd],axis=1)
		n += 1

F1_median.to_csv("/mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/AucROC_PCC/Naive_median_F1_test.txt",index=True, header=True,sep="\t")
F1_max.to_csv("/mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/AucROC_PCC/Naive_max_F1_test.txt",index=True, header=True,sep="\t")


