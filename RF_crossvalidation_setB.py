import sys, os
import pandas as pd
import numpy as np
from datetime import datetime
import time
import ML_function_for_CV as ML
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import cross_val_predict

def Performance_MC(y, pred, classes):
	from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
	cm = confusion_matrix(y, pred, labels=classes)
	accuracy = accuracy_score(y, pred)
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
	macro_f1 = np.mean(f1)
	return {'cm':cm, 'accuracy':accuracy,'macro_f1':macro_f1,'f1_MC':f1}

def DefineClf_RandomForest(n_estimators,max_depth,max_features):
	from sklearn.ensemble import RandomForestClassifier
	clf = RandomForestClassifier(n_estimators=int(n_estimators), 
		max_depth=max_depth, 
		max_features=max_features,
		criterion='gini')
	return clf

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
		
	# grid search 	
	CV_performance = {}	# CV_performance [max_depth,max_features,n_estimators]=[F1_macro,F1_macro,...]
	parameters = {'max_depth':[3, 5, 10], 'max_features': [0.1, 0.25, 0.5, 0.75, 'sqrt', 'log2', None], 'n_estimators': [100,500,1000]}
	for max_depth in parameters['max_depth']:
		for max_features in parameters['max_features']:
			for n_estimators in parameters['n_estimators']:
				clf = DefineClf_RandomForest(n_estimators,max_depth,max_features)
				for cv_number in range(1,6):
					df = pd.read_csv(path+DF+ '_CV_%s_features.txt'%cv_number, sep='\t', index_col = 0)
					df_test = df[df.index.isin(test)]
					y_test = df_test['Class']
					X_test = df_test.drop(['Class'], axis=1)
					test_classes = y_test.unique()
					#split genes into traning, validiating, testing subset
					with open('Genes_for_5_training_set%s.txt'%cv_number) as train_file:
						train = train_file.read().splitlines()
					with open('Genes_for_5_validation_set%s.txt'%cv_number) as validation_file:
						validation = validation_file.read().splitlines()
					df_train = df[df.index.isin(train)]
					df_validation = df[df.index.isin(validation)]
					y_train = df_train['Class']
					X_train = df_train.drop(['Class'], axis=1)
					y_validation = df_validation['Class']
					X_validation = df_validation.drop(['Class'], axis=1)
					classes = y_validation.unique()
					clf.fit(X_train,y_train)
					cv_proba = clf.predict_proba(df_validation.drop(['Class'], axis=1))
					cv_pred = clf.predict(df_validation.drop(['Class'], axis=1))
					result = Performance_MC(y_validation, cv_pred, classes)
					if 'macro_f1' in result:
						if '%s_%s_%s'%(n_estimators, max_depth, max_features) not in CV_performance:
							CV_performance['%s_%s_%s'%(n_estimators, max_depth, max_features)] = []				
						CV_performance['%s_%s_%s'%(n_estimators, max_depth, max_features)].append(result['macro_f1']) 

	# get the average F1 for each parameter set
	CV_p = {}
	for parameter in CV_performance:
		if np.mean(CV_performance[parameter]) not in CV_p:
			CV_p[np.mean(CV_performance[parameter])] = []
		CV_p[np.mean(CV_performance[parameter])].append(parameter)
		
	parameter = CV_p[max(CV_p.keys())][0]
	n_estimators,max_depth,max_features = parameter.split('_')
	n_estimators = int(n_estimators)
	max_depth = int(max_depth)
	print(max_features)
	if max_features not in ['sqrt', 'log2', None,'auto']:
		try:
			max_features = float(max_features)
		except:
			print('max_features is not a float value')
	if max_features == 'None':
		max_features = None
	
	####### Run ML models #######
	## Make empty dataframes
	conf_matrices = pd.DataFrame(columns = np.insert(arr = classes.astype(np.str), obj = 0, values = 'Class'))
	imp = pd.DataFrame(index = list(df.drop(['Class'], axis=1)))
	accuracies = []
	f1_array = np.array([np.insert(arr = classes.astype(np.str), obj = 0, values = 'M')])
	accuracies_ho = []
	f1_array_ho = np.array([np.insert(arr = test_classes.astype(np.str), obj = 0, values = 'M')])
	results = []
	results_ho = []
	df_all = df.copy()
	df_proba = pd.DataFrame(data=df_all['Class'], index=df_all.index, columns=['Class'])
	clf = DefineClf_RandomForest(n_estimators,max_depth,max_features)
	for cv_number in range(1,6):
		df = pd.read_csv(path+DF+ '_CV_%s_features.txt'%cv_number, sep='\t', index_col = 0)
		df_test = df[df.index.isin(test)]
		y_test = df_test['Class']
		X_test = df_test.drop(['Class'], axis=1)
		#split genes into traning, validiating, testing subset
		with open('Genes_for_5_training_set%s.txt'%cv_number) as train_file:
			train = train_file.read().splitlines()
		with open('Genes_for_5_validation_set%s.txt'%cv_number) as validation_file:
			validation = validation_file.read().splitlines()
		df_train = df[df.index.isin(train)]
		df_validation = df[df.index.isin(validation)]
		y_train = df_train['Class']
		X_train = df_train.drop(['Class'], axis=1)
		y_validation = df_validation['Class']
		X_validation = df_validation.drop(['Class'], axis=1)
		classes = y_validation.unique()
		test_classes = y_test.unique()
		clf.fit(X_train,y_train)
		importances = clf.feature_importances_
		imp[cv_number] = importances
		cv_proba = clf.predict_proba(df_validation.drop(['Class'], axis=1))
		cv_pred = clf.predict(df_validation.drop(['Class'], axis=1))
		result = Performance_MC(y_validation, cv_pred, classes)
		ho_proba = clf.predict_proba(df_test.drop(['Class'], axis=1))
		ho_pred = clf.predict(df_test.drop(['Class'], axis=1))
		ho_result = Performance_MC(y_test, ho_pred, test_classes)
		score_columns = []
		for clss in classes:
			score_columns.append("%s_score_%s"%(clss,cv_number))
		df_sel_scores = pd.DataFrame(data=cv_proba,index=df_validation.index,columns=score_columns)
		df_ho_scores = pd.DataFrame(data=ho_proba,index=df_test.index,columns=score_columns)
		current_scores =  pd.concat([df_sel_scores,df_ho_scores], axis = 0)
		if 'accuracy' in result:
			accuracies.append(result['accuracy'])
		if 'macro_f1' in result:
			f1_temp_array = np.insert(arr = result['f1_MC'], obj = 0, values = result['macro_f1'])
			f1_array = np.append(f1_array, [f1_temp_array], axis=0)
		if 'cm' in result:
			cmatrix = pd.DataFrame(result['cm'], columns = classes)
			cmatrix.insert(0,"Class",classes,True)
			conf_matrices = pd.concat([conf_matrices, cmatrix],axis=0)
		if 'accuracy' in ho_result:
			accuracies_ho.append(ho_result['accuracy'])
		if 'macro_f1' in ho_result:
			ho_f1_temp_array = np.insert(arr = ho_result['f1_MC'], obj = 0, values = ho_result['macro_f1'])
			f1_array_ho = np.append(f1_array_ho, [ho_f1_temp_array], axis=0)
		results.append(result)
		df_proba = pd.concat([df_proba,current_scores], axis = 1,sort=True)

###### save the Output ######
	if len(classes) > 2:
		summary_cols = []
		mc_score_columns = []
		keep_for_summary = ['Class', 'Prediction']
		for class_nm in reversed(classes): # get std
			class_proba_cols = [c for c in df_proba.columns if c.startswith(class_nm+'_score_')]
			df_proba.insert(loc=1, column = class_nm+'_score_stdev', value = df_proba[class_proba_cols].std(axis=1))
			summary_cols.insert(0,class_nm+'_score_stdev')

		for class_nm in reversed(classes): # get median
			summary_cols.insert(0,class_nm+'_score_Median')
			mc_score_columns.append(class_nm+'_score_Median')
			keep_for_summary.append(class_nm + '_score_Median')
			class_proba_cols = [c for c in df_proba.columns if c.startswith(class_nm+'_score_')]
			df_proba.insert(loc=1, column = class_nm+'_score_Median', value = df_proba[class_proba_cols].median(axis=1))
			
		
		# Find the max mc_score and set to Prediction column (remove the _score_Median string)
		df_proba.insert(loc=1, column = 'Prediction', value = df_proba[mc_score_columns].idxmax(axis=1))
		df_proba['Prediction'] = df_proba.Prediction.str.replace('_score_Median','')

		
		# Count the number of times an instance of class x is predicted as class y 		
		summary_df_proba = df_proba[['Class', 'Prediction',class_nm + '_score_Median']].groupby(['Class', 'Prediction']).agg('count').unstack(level=1)
		summary_df_proba.columns = summary_df_proba.columns.droplevel()

		# Check to make sure each class was predicted at least once
		for cl in classes:
			if cl not in list(summary_df_proba):
				print('No instances were predicted as class: %s' % cl)
				summary_df_proba[cl] = 0
		summary_df_proba['n_total'] = summary_df_proba[classes].sum(axis=1)


		for class_nm in classes:
			summary_df_proba[str(class_nm) + '_perc'] = summary_df_proba[class_nm]/summary_df_proba['n_total']

		# save scores file
		scores_file = save_path + DF + "_RF_scores.txt"
		out_scores = open(scores_file,"w")
		out_scores.write("#ID\t"+pd.DataFrame.to_csv(df_proba[["Class"]+summary_cols],sep="\t").strip()+"\n")
		out_scores.write("#ID\t"+pd.DataFrame.to_csv(df_proba,sep="\t").strip()+"\n")
		out_scores.close()
		
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

		print("\nML Results: \nAccuracy: %03f (+/- stdev %03f)\nF1 (macro): %03f (+/- stdev %03f)\n" % (
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
		print("\nML Results from Hold Out Validation: \nAccuracy: %03f (+/- stdev %03f)\nF1 (macro): %03f (+/- stdev %03f)\n" % (AC_ho, AC_std_ho, MacF1_ho, MacF1_std_ho))

		# Save detailed results file 
		n_features = df.shape[1] - 1 
		out = open(save_path + DF + "_RF_results.txt", 'w')
		out.write('ID: %s\nAlgorithm: RF\nTrained on classes: %s\nNumber of features: %i\n' % (DF, classes, n_features))
		out.write('Parameters used:n_estimators:%s\nmax_depth:%s\nmax_features:%s\n' % (n_estimators,max_depth,max_features))

		out.write('\nMetric\tMean\tSD\nAccuracy\t%05f\t%05f\nF1_macro\t%05f\t%05f\n' % (AC, AC_std, MacF1, MacF1_std))
		for cla in f1.columns:
			if 'M_F1' not in cla:
				out.write('%s\t%05f\t%05f\n' % (cla, np.mean(f1[cla]), np.std(f1[cla])))
		
		out.write('\nMean Balanced Confusion Matrix:\n')
		cm_mean.to_csv(out, mode='a', sep='\t')
		out.write('\n\nCount and percent of instances of each class (row) predicted as a class (col):\n')
		summary_df_proba.to_csv(out, mode='a', header=True, sep='\t')
			
		# Add results from hold out
		out.write('\n\nResults from the hold out validation set\n')
		out.write('HO Accuracy\t%05f +/-%05f\nHO F1_macro\t%05f +/-%05f\n' % (AC_ho, AC_std_ho, MacF1_ho, MacF1_std_ho))
		for cla in f1_ho.columns:
			if 'M_F1' not in cla:
				out.write('%s\t%05f\t%05f\n' % (cla, np.mean(f1_ho[cla]), np.std(f1_ho[cla])))

		out.close()

		imp['mean_imp'] = imp.mean(axis=1)
		imp = imp.sort_values('mean_imp', 0, ascending = False)
		imp_out = save_path + DF + "_imp"
		imp['mean_imp'].to_csv(imp_out, sep = "\t", index=True)
	

if __name__ == '__main__':
	main()
