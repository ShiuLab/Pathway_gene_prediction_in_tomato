import sys,os,argparse
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from imblearn.over_sampling import SMOTE
import joblib

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn


def Performance_MC(y, pred, classes):
	from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
	cm = confusion_matrix(y, pred, labels=classes)
	accuracy = accuracy_score(y, pred)
	pred = pd.DataFrame(pred,index=y.index)
	df_tem = pd.concat([y,pred],axis=1)
	df_tem.columns = ['y_true','y_pred']
	f1 = []
	for y_list in classes:
		P = df_tem[df_tem.y_true==y_list]
		TP = P[P.y_pred == y_list]
		FN = P.shape[0] - TP.shape[0]
		N = df_tem[df_tem.y_true!=y_list]
		FP = N[N.y_pred == y_list]
		f1_tem = Fmeasure(TP.shape[0],FP.shape[0],FN)
		f1.append(f1_tem)
	macro_f1 = np.mean(f1)
	return {'cm':cm, 'accuracy':accuracy,'macro_f1':macro_f1,'f1_MC':f1}

def DefineClf_AdaBoost(base_estimator,learning_rate):
	clf = AdaBoostClassifier(base_estimator, 
		learning_rate=learning_rate,
		n_estimators=50, 
		random_state=42)
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
	parser = argparse.ArgumentParser(description='This code is for building the AdaBoostClassifier models, where the clf_estimator is the RF model with the best hyperparameters for each dataset. Numbers of genes in pathways are balanced for the training set')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-df_short_name', help='feature matrix, for Set B, use the short name, for Set A, use the full name of the expression matrix', required=True)
	req_group.add_argument('-path', help='path to the feature matrix', required=True)
	req_group.add_argument('-save_path', help='path to save the outputs', required=True)
	req_group.add_argument('-test_gene_list', help='Genes_for_testing.txt', required=True)
	req_group.add_argument('-train_gene_list', help='Genes_for_training.txt', required=True)
	req_group.add_argument('-dataset', help='setA or setB', required=True)
	req_group.add_argument('-clf_estimator', help='', required=True)
	req_group.add_argument('-max_depth', help='', required=True)
	req_group.add_argument('-max_features', help='', required=True)
	req_group.add_argument('-n_estimators', help='', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	

	DF = args.df_short_name
	path = args.path
	save_path = args.save_path
	TEST = args.test_gene_list
	TRAIN = args.train_gene_list
	dataset = args.dataset
	clf_estimator = args.clf_estimator
	max_depth = int(args.max_depth)
	max_features = args.max_features
	n_estimators = int(args.n_estimators)

	if max_features in ['0.1', '0.25', '0.5', '0.75']:
		max_features = float(max_features)
	elif max_features == 'None':
		max_features = None
		
	with open(TEST) as test_file:
		test = test_file.read().splitlines()

	with open(TRAIN) as training_file:
		training = training_file.read().splitlines()
		
	if dataset == 'setB':
		df = pd.read_csv(path+DF+ '_CV_1_features.txt', sep='\t', index_col = 0)
		short_name = DF
	if dataset == 'setA':
		expression = pd.read_csv(path+DF, sep='\t', index_col = 0)
		pathway_annotation = pd.read_csv('Sly_pathway_annotation_20190117_with_expression_5_members_nonoverlapping.txt', sep='\t', index_col = 1,header=None)
		pathway_annotation.columns = ['Class']
		df = pd.concat([pathway_annotation,expression],axis = 1)
		short_name = DF
		
	df_test = df[df.index.isin(test)]
	y_test = df_test['Class']
	X_test = df_test.drop(['Class'], axis=1)
	df_training = df[df.index.isin(training)]
	y_training = df_training['Class']
	X_training = df_training.drop(['Class'], axis=1)
	# grid search 
	y = df['Class']
	classes = y.unique()
	df_test = df[df.index.isin(test)]
	y_test = df_test['Class']
	test_classes = y_test.unique()
	CV_performance = {}	# CV_performance [max_depth,max_features,n_estimators]=[F1_macro,F1_macro,...]
	out = open(save_path + short_name + '_%s_CV_performance.txt'%clf_estimator,'w')
	perf = []
	if clf_estimator == 'RF':
		base_estimator = RandomForestClassifier(n_estimators=int(n_estimators), 
		max_depth=max_depth, 
		max_features=max_features,
		criterion='gini',
		random_state=42)
	for learning_rate in [10,1,0.01,0.001]:
		clf = DefineClf_AdaBoost(base_estimator,learning_rate)
		for cv_number in range(1,6):
			if dataset == 'setB':
				df = pd.read_csv(path+DF+ '_CV_%s_features.txt'%cv_number, sep='\t', index_col = 0)
			#split genes into traning, validiating, testing subset
			with open('Genes_for_5_training_set%s.txt'%cv_number) as train_file:
				train = train_file.read().splitlines()
			with open('Genes_for_5_validation_set%s.txt'%cv_number) as validation_file:
				validation = validation_file.read().splitlines()
			df_train = df[df.index.isin(train)]
			df_validation = df[df.index.isin(validation)]
			y_train = df_train['Class']
			X_train = df_train.drop(['Class'], axis=1)
			# upsample the minorities
			smote = SMOTE(sampling_strategy='not majority',random_state=42,k_neighbors=3)
			X_train, y_train = smote.fit_sample(X_train, y_train)
			y_validation = df_validation['Class']
			X_validation = df_validation.drop(['Class'], axis=1)
			clf.fit(X_train,y_train)
			cv_pred = clf.predict(df_validation.drop(['Class'], axis=1))
			result = Performance_MC(y_validation, cv_pred, classes)
			if 'macro_f1' in result:
				if '%s__%s__%s__%s___%s'%(clf_estimator,n_estimators,max_depth,max_features,learning_rate) not in CV_performance:
					CV_performance['%s__%s__%s__%s___%s'%(clf_estimator,n_estimators,max_depth,max_features,learning_rate)] = []				
				CV_performance['%s__%s__%s__%s___%s'%(clf_estimator,n_estimators,max_depth,max_features,learning_rate)].append(result['macro_f1']) 

	# get the average F1 for each parameter set
	CV_p = {}
	for parameter in CV_performance:
		if np.mean(CV_performance[parameter]) not in CV_p:
			CV_p[np.mean(CV_performance[parameter])] = []
		CV_p[np.mean(CV_performance[parameter])].append(parameter)

	parameter = CV_p[max(CV_p.keys())][0]
	clf_estimator,n_estimators,max_depth,max_features = parameter.split('___')[0].split('__')
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
	learning_rate = parameter.split('___')[1]
	if learning_rate in ['10','1']:
		learning_rate = int(learning_rate)
	else:
		learning_rate = float(learning_rate)

	
	####### Run ML models #######
	## Make empty dataframes
	conf_matrices = pd.DataFrame(columns = np.insert(arr = classes.astype(np.str), obj = 0, values = 'Class'))
	# imp = pd.DataFrame(index = list(df.drop(['Class'], axis=1)))
	accuracies = []
	accuracies_ho = []
	f1_array = np.array([np.insert(arr = classes.astype(np.str), obj = 0, values = 'M')])
	accuracies_ho = []
	accuracies_ho = []
	f1_array_ho = np.array([np.insert(arr = test_classes.astype(np.str), obj = 0, values = 'M')])
	results = []
	results_ho = []
	df_all = df.copy()
	df_proba_test = pd.DataFrame(data=df_test['Class'], index=df_test.index, columns=['Class'])
	if clf_estimator == 'RF':
		base_estimator = RandomForestClassifier(n_estimators=int(n_estimators), 
		max_depth=max_depth, 
		max_features=max_features,
		criterion='gini',
		random_state=42)
	clf = DefineClf_AdaBoost(base_estimator,learning_rate)
	Prediction_testing = pd.DataFrame(y_test)
	for cv_number in range(1,6):
		if dataset == 'setB':
			df = pd.read_csv(path+DF+ '_CV_%s_features.txt'%cv_number, sep='\t', index_col = 0)
		#split genes into traning, validiating, testing subset
		with open('Genes_for_5_training_set%s.txt'%cv_number) as train_file:
			train = train_file.read().splitlines()
		with open('Genes_for_5_validation_set%s.txt'%cv_number) as validation_file:
			validation = validation_file.read().splitlines()
		df_train = df[df.index.isin(train)]
		df_validation = df[df.index.isin(validation)]
		y_train = df_train['Class']
		X_train = df_train.drop(['Class'], axis=1)
		# upsample the minorities
		smote = SMOTE(sampling_strategy='not majority',random_state=42,k_neighbors=3)
		X_train, y_train = smote.fit_sample(X_train, y_train)
		df_training = pd.concat([y_train,X_train],axis=1)
		y_validation = df_validation['Class']
		X_validation = df_validation.drop(['Class'], axis=1)
		clf.fit(X_train,y_train)
		# importances = clf.feature_importances_
		# imp[cv_number] = importances
		
		proba = clf.predict_proba(df_validation.drop(['Class'], axis=1))
		pred = clf.predict(df_validation.drop(['Class'], axis=1))
		joblib.dump(clf, save_path + short_name + "_AdaBoost_SMOTE_%s_%s.pkl"%(cv_number,dataset))

		pred_va = pd.DataFrame(y_validation)
		pred_va['Prediction'] = pred
		if cv_number == 1:
			Prediction_validation = pred_va
		else:
			Prediction_validation = pd.concat([Prediction_validation,pred_va],axis=0)
		result = Performance_MC(y_validation, pred, classes)
		
		ho_proba = clf.predict_proba(df_test.drop(['Class'], axis=1))
		ho_pred = clf.predict(df_test.drop(['Class'], axis=1))
		Prediction_testing['Prediction_%s'%cv_number] = ho_pred
		results_ho = Performance_MC(y_test, ho_pred, test_classes)
		
		score_columns = []
		for clss in classes:
			score_columns.append("%s_score_%s"%(clss,cv_number))
		df_sel_scores = pd.DataFrame(data=proba,index=df_validation.index,columns=score_columns)
		df_proba_col = pd.DataFrame(data=df_validation['Class'], index=df_validation.index, columns=['Class'])
		df_proba_tem = pd.concat([df_proba_col,df_sel_scores], axis = 1, sort = True)
		if cv_number==1:
			df_proba = df_proba_tem
		else:
			df_proba = pd.concat([df_proba,df_proba_tem], axis = 0)
		df_ho_scores = pd.DataFrame(data=ho_proba,index=df_test.index,columns=score_columns)
		if 'accuracy' in result:
			accuracies.append(result['accuracy'])
		if 'macro_f1' in result:
			f1_temp_array = np.insert(arr = result['f1_MC'], obj = 0, values = result['macro_f1'])
			f1_array = np.append(f1_array, [f1_temp_array], axis=0)
		if 'cm' in result:
			cmatrix = pd.DataFrame(result['cm'], columns = classes)
			cmatrix.insert(0,"Class",classes,True)
			conf_matrices = pd.concat([conf_matrices, cmatrix],axis=0)
		if 'accuracy' in results_ho:
			accuracies_ho.append(results_ho['accuracy'])
		if 'macro_f1' in results_ho:
			ho_f1_temp_array = np.insert(arr = results_ho['f1_MC'], obj = 0, values = results_ho['macro_f1'])
			f1_array_ho = np.append(f1_array_ho, [ho_f1_temp_array], axis=0)
		df_proba_test = pd.concat([df_proba_test,df_ho_scores], axis = 1,sort=True)

###### save the Output ######
	if len(classes) > 2:
		summary_cols = []
		mc_score_columns = []
		keep_for_summary = ['Class', 'Prediction']
		for class_nm in reversed(classes): # get std
			class_proba_cols = [c for c in df_proba.columns if c.startswith(class_nm+'_score_')]
			df_proba.insert(loc=1, column = class_nm+'_score_stdev', value = df_proba[class_proba_cols].std(axis=1))
			
			class_proba_cols_test = [c for c in df_proba_test.columns if c.startswith(class_nm+'_score_')]
			df_proba_test.insert(loc=1, column = class_nm+'_score_stdev', value = df_proba_test[class_proba_cols_test].std(axis=1))
			summary_cols.insert(0,class_nm+'_score_stdev')

		for class_nm in reversed(classes): # get median
			summary_cols.insert(0,class_nm+'_score_Median')
			mc_score_columns.append(class_nm+'_score_Median')
			keep_for_summary.append(class_nm + '_score_Median')
			class_proba_cols = [c for c in df_proba.columns if c.startswith(class_nm+'_score_')]
			df_proba.insert(loc=1, column = class_nm+'_score_Median', value = df_proba[class_proba_cols].median(axis=1))
			
			class_proba_cols_test = [c for c in df_proba_test.columns if c.startswith(class_nm+'_score_')]
			df_proba_test.insert(loc=1, column = class_nm+'_score_Median', value = df_proba_test[class_proba_cols_test].median(axis=1))
			
		
		# Find the max mc_score and set to Prediction column (remove the _score_Median string)
		df_proba.insert(loc=1, column = 'Prediction', value = df_proba[mc_score_columns].idxmax(axis=1))
		df_proba['Prediction'] = df_proba.Prediction.str.replace('_score_Median','')

		df_proba_test.insert(loc=1, column = 'Prediction', value = df_proba_test[mc_score_columns].idxmax(axis=1))
		df_proba_test['Prediction'] = df_proba_test.Prediction.str.replace('_score_Median','')

		
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
		scores_file = save_path + short_name + "_AdaBoost_%s_validation_scores.txt"%dataset
		out_scores = open(scores_file,"w")
		out_scores.write("#ID\t"+pd.DataFrame.to_csv(df_proba,sep="\t").strip()+"\n")
		out_scores.close()
		
		scores_file = save_path + short_name + "_AdaBoost_%s_test_scores.txt"%dataset
		out_scores = open(scores_file,"w")
		out_scores.write("#ID\t"+pd.DataFrame.to_csv(df_proba_test,sep="\t").strip()+"\n")
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
		out = open(save_path + short_name + "_AdaBoost_%s_results.txt"%dataset, 'w')
		out.write('ID: %s\nAlgorithm: NN\nTrained on classes: %s\nNumber of features: %i\n' % (DF, classes, n_features))
		out.write('Parameters used:clf_estimator:%s\nmax_depth:%s\nmax_features:%s\nn_estimators:%s\nlearning_rate:%s\n' % (clf_estimator,max_depth,max_features,n_estimators,learning_rate))

		out.write('\n\nResults for prediction on validation set:\n')
		out.write('\nMetric\tMean\tSD\nAccuracy\t%05f\t%05f\nF1_macro\t%05f\t%05f\n' % (AC, AC_std, MacF1, MacF1_std))
		for cla in f1.columns:
			if 'M_F1' not in cla:
				out.write('%s\t%05f\t%05f\n' % (cla, np.mean(f1[cla]), np.std(f1[cla])))
		
		# Add results from hold out
		out.write('\n\nResults from the hold out set:\n')
		out.write('HO Accuracy\t%05f +/-%05f\nHO F1_macro\t%05f +/-%05f\n' % (AC_ho, AC_std_ho, MacF1_ho, MacF1_std_ho))
		for cla in f1_ho.columns:
			if 'M_F1' not in cla:
				out.write('%s\t%05f\t%05f\n' % (cla, np.mean(f1_ho[cla]), np.std(f1_ho[cla])))

		out.write('\n\nMean Balanced Confusion Matrix:\n')
		cm_mean.to_csv(out, mode='a', sep='\t')
		out.write('\n\nCount and percent of instances of each class (row) predicted as a class (col):\n')
		summary_df_proba.to_csv(out, mode='a', header=True, sep='\t')

		out.close()

		# imp['mean_imp'] = imp.mean(axis=1)
		# imp = imp.sort_values('mean_imp', 0, ascending = False)
		# imp_out = save_path + short_name + "_%s_imp"%dataset
		# imp['mean_imp'].to_csv(imp_out, sep = "\t", index=True)
		
		Prediction_validation.to_csv(save_path + short_name + "_AdaBoost_%s_validation_prediction.txt"%dataset,index=True, header=True,sep="\t")
		Prediction_testing.to_csv(save_path + short_name + "_AdaBoost_%s_testing_prediction.txt"%dataset,index=True, header=True,sep="\t")
		

	

if __name__ == '__main__':
	main()
