import sklearn
import pandas as pd
import numpy as np
import sys,os
from sklearn.cluster import KMeans
import scipy
import scipy.stats as stats
import math

def Performance_MC(y, pred, classes):
	from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
	cm = confusion_matrix(y, pred, labels=classes)
	accuracy = accuracy_score(y, pred)
	df_tem = pd.DataFrame([y,pred]).T
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
	macro_f1 = np.mean(f1)
	return {'cm':cm, 'accuracy':accuracy,'macro_f1':macro_f1,'f1_MC':f1}

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


def Enrichment_clustering(cluster_result,n_clusters):
	Enrichment_C_P = {}
	Enrichment = {}
	for pathway in np.unique(cluster_result.Class):
		for cluster in range(0,n_clusters):
			P = cluster_result[cluster_result.Class == pathway]
			C_P = P[P.Cluster == cluster]
			C = cluster_result[cluster_result.Cluster == cluster]
			if C_P.shape[0] != 0:
				logRatio =  math.log((float(C_P.shape[0])/C.shape[0])/(float(P.shape[0])/cluster_result.shape[0]))
				pvalue = stats.fisher_exact([[C_P.shape[0], P.shape[0] - C_P.shape[0]], [C.shape[0] - C_P.shape[0], cluster_result.shape[0] - C.shape[0] - P.shape[0] + C_P.shape[0]]])[1]
				if pvalue < 0.05:
					if cluster not in Enrichment:
						Enrichment[cluster] = {}
					if logRatio not in Enrichment[cluster]:
						Enrichment[cluster][logRatio] = []
						Enrichment[cluster][logRatio].append([pathway,pvalue])
			for cluster in Enrichment:
				best_pa_tem = Enrichment[cluster][max(Enrichment[cluster].keys())]
				if len(best_pa_tem) == 1:
					best_pa = best_pa_tem[0][0]
				else:
					min_p = 1
					for enrichment in best_pa_tem:
						min_p = min(min_p,enrichment[1])
					for enrichment in best_pa_tem:
						if enrichment[1] == min_p:
							best_pa = enrichment[0]
				Enrichment_C_P[cluster] = best_pa
	return(Enrichment_C_P)

def main():
	for i in range (1,len(sys.argv),2):
		if sys.argv[i].lower() == "-df_short_name":
			DF = sys.argv[i+1]
		if sys.argv[i].lower() == "-path": ### path to feature files
			path = sys.argv[i+1]
		if sys.argv[i].lower() == "-save_path": ### path to feature files
			save_path = sys.argv[i+1]
	
	TEST = 'Genes_for_testing.txt'
	with open(TEST) as test_file:
		test = test_file.read().splitlines()
	
	expression = pd.read_csv(path+DF, sep='\t', index_col = 0)
	pathway_annotation = pd.read_csv('Sly_pathway_annotation_20190117_with_expression_5_members_nonoverlapping.txt', sep='\t', index_col = 1,header=None)
	pathway_annotation.columns = ['Class']
	df = pd.concat([pathway_annotation,expression],axis = 1)
	print(df.columns)
	y = df['Class']
	classes = y.unique()
	df_test = df[df.index.isin(test)]
	y_test = df_test['Class']
	test_classes = y_test.unique()
	conf_matrices = pd.DataFrame(columns = np.insert(arr = classes.astype(np.str), obj = 0, values = 'Class'))
	accuracies = []
	f1_array = np.array([np.insert(arr = classes.astype(np.str), obj = 0, values = 'M')])
	accuracies_ho = []
	f1_array_ho = np.array([np.insert(arr = test_classes.astype(np.str), obj = 0, values = 'M')])
	results = []
	df_all = df.copy()
	df_proba = pd.DataFrame(data=df_all['Class'], index=df_all.index, columns=['Class'])
	for n_clusters in [5,10,25,50,85,100,200,300,400,500]:
		for cv_number in range(1,6):
			df_test = df[df.index.isin(test)]
			y_test = df_test['Class']
			test_classes = y_test.unique()
			X_test = df_test.drop(['Class'], axis=1)
			with open('Genes_for_5_training_set%s.txt'%cv_number) as train_file:
				train = train_file.read().splitlines()
			with open('Genes_for_5_validation_set%s.txt'%cv_number) as validation_file:
				validation = validation_file.read().splitlines()
			df_train = df[df.index.isin(train)]
			df_validation = df[df.index.isin(validation)]
			X_train = df_train.drop(['Class'], axis=1)
			X_validation = df_validation.drop(['Class'], axis=1)
			y_train = df_train['Class']
			y_validation = df_validation['Class']
			mat = X_train.as_matrix()  # Convert DataFrame to matrix
			mat_validation = X_validation.as_matrix()
			mat_test = X_test.as_matrix()
			km = sklearn.cluster.KMeans(n_clusters=n_clusters,n_init=3, n_jobs=5,max_iter=500) # Using sklearn
			km.fit(mat)
			train_labels = km.labels_   # Get cluster assignment labels
			train_tem = pd.DataFrame([train_labels]).T # Format results as a DataFrame
			train_tem.index = X_train.index
			train_tem.columns = ['Cluster']
			train_res = pd.concat([y_train,train_tem],axis=1)
			E_C_P = Enrichment_clustering(train_res,n_clusters)
			cv_labels = km.predict(mat_validation)
			cv_tem = pd.DataFrame([cv_labels]).T
			cv_tem.index = X_validation.index
			cv_tem.columns = ['Cluster']
			cv_res = pd.concat([y_validation,cv_tem],axis=1)
			for i in range(0,cv_res.shape[0]):
				try:
					cv_res.iloc[i,1] = E_C_P[cv_res.iloc[i,1]]
				except:
					cv_res.iloc[i,1] = '%s'%cv_res.iloc[i,1]
					print('%s was not enriched for any pathway'%cv_res.iloc[i,1])
			result = Performance_MC(cv_res.Class, cv_res.Cluster, classes)
			if 'accuracy' in result:
				accuracies.append(result['accuracy'])
			if 'macro_f1' in result:
				f1_temp_array = np.insert(arr = result['f1_MC'], obj = 0, values = result['macro_f1'])
				f1_array = np.append(f1_array, [f1_temp_array], axis=0)
			test_labels = km.predict(mat_test)
			test_tem = pd.DataFrame([test_labels]).T
			test_tem.index = X_test.index
			test_tem.columns = ['Cluster']
			test_res = pd.concat([y_test,test_tem],axis=1)
			for i in range(0,test_res.shape[0]):
				try:
					test_res.iloc[i,1] = E_C_P[test_res.iloc[i,1]]
				except:
					test_res.iloc[i,1] = '%s'%test_res.iloc[i,1]
					print('%s was not enriched for any pathway'%test_res.iloc[i,1])
			ho_result = Performance_MC(test_res.Class, test_res.Cluster, test_classes)
			if 'accuracy' in ho_result:
				accuracies_ho.append(ho_result['accuracy'])
			if 'macro_f1' in ho_result:
				ho_f1_temp_array = np.insert(arr = ho_result['f1_MC'], obj = 0, values = ho_result['macro_f1'])
				f1_array_ho = np.append(f1_array_ho, [ho_f1_temp_array], axis=0)
			results.append(result)

		f1 = pd.DataFrame(f1_array)
		f1.columns = f1.iloc[0]
		f1 = f1[1:]
		f1.columns = [str(col) + '_F1' for col in f1.columns]
		f1 = f1.astype(float)		
		
		# Calculate accuracy and f1 stats
		AC = np.mean(accuracies)
		AC_std = np.std(accuracies)
		MacF1 = f1['M_F1'].mean()
		MacF1_std = f1['M_F1'].std()

		print("\nCluster Results: \nAccuracy: %03f (+/- stdev %03f)\nF1 (macro): %03f (+/- stdev %03f)\n" % (
		AC, AC_std, MacF1, MacF1_std))
		
		for class_name in classes:
			conf_matrices[class_name] = conf_matrices[class_name].astype(float)
		cm_mean = conf_matrices.groupby('Class').mean()

		# Unpack results from hold out
		f1_ho = pd.DataFrame(f1_array_ho)
		f1_ho.columns = f1_ho.iloc[0]
		f1_ho = f1_ho[1:]
		f1_ho.columns = [str(col) + '_F1' for col in f1_ho.columns]
		f1_ho = f1_ho.astype(float)	
		AC_ho = np.mean(accuracies_ho)
		AC_std_ho = np.std(accuracies_ho)
		MacF1_ho = f1_ho['M_F1'].mean()
		MacF1_std_ho = f1_ho['M_F1'].std()
		print("\nCluster Results from Hold Out Validation: \nAccuracy: %03f (+/- stdev %03f)\nF1 (macro): %03f (+/- stdev %03f)\n" % (AC_ho, AC_std_ho, MacF1_ho, MacF1_std_ho))

		# Save detailed results file 
		n_features = df.shape[1] - 1 
		out = open(save_path+DF + "_Kmean_setA_%s_results.txt"%n_clusters, 'w')

		out.write('\nMetric\tMean\tSD\nAccuracy\t%05f\t%05f\nF1_macro\t%05f\t%05f\n' % (AC, AC_std, MacF1, MacF1_std))
		for cla in f1.columns:
			if 'M_F1' not in cla:
				out.write('%s\t%05f\t%05f\n' % (cla, np.mean(f1[cla]), np.std(f1[cla])))
		
		
		# Add results from hold out
		out.write('\n\nResults from the hold out validation set\n')
		out.write('HO Accuracy\t%05f +/-%05f\nHO F1_macro\t%05f +/-%05f\n' % (AC_ho, AC_std_ho, MacF1_ho, MacF1_std_ho))
		for cla in f1_ho.columns:
			if 'M_F1' not in cla:
				out.write('%s\t%05f\t%05f\n' % (cla, np.mean(f1_ho[cla]), np.std(f1_ho[cla])))

		out.close()


if __name__ == '__main__':
	main()

