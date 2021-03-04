# source /mnt/home/mengfanr/run_local_tensorflow_module_load_gpu_2.2.sh
import sys,os,argparse
import pandas as pd
import numpy as np
from imblearn.over_sampling import SMOTE
import tensorflow as tf
from tensorflow import keras
from sklearn.metrics import f1_score
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
	parser = argparse.ArgumentParser(description='This code is for building the keras Sequential API based model, where numbers of genes in pathways are balanced for the training set.')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-df_short_name', help='feature matrix, for Set B, use the short name, for Set A, use the full name of the expression matrix', required=True)
	req_group.add_argument('-path', help='path to the feature matrix', required=True)
	req_group.add_argument('-save_path', help='path to save the outputs', required=True)
	req_group.add_argument('-test_gene_list', help='Genes_for_testing.txt', required=True)
	req_group.add_argument('-train_gene_list', help='Genes_for_training.txt', required=True)
	req_group.add_argument('-dataset', help='setA or setB', required=True)

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
		
	np.random.seed(42)
	tf.random.set_seed(42)
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
	parameters = {'activation':['relu','sigmoid','selu','elu'], 'kernel_initializer': ['RandomNormal','RandomUniform','VarianceScaling','Zeros','GlorotNormal','Orthogonal'], 'Dropout_rate': [0.1,0.2,0.5]}
	for activation in parameters['activation']:
		for kernel_initializer in parameters['kernel_initializer']:
			for Dropout_rate in parameters['Dropout_rate']:
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
					### need to convert the string classes to numbers
					if cv_number == 1:
						Name = {}
						Rev = {}
						n = 0
						for i in range(0,y_train.shape[0]):
							if y_train.iloc[i] not in Name:
								Name[y_train.iloc[i]] = n
								Rev[n] = y_train.iloc[i]
								y_train.iloc[i] = n
								n += 1
							else:
								y_train.iloc[i] = Name[y_train.iloc[i]]
						f = open("dict0.txt","w")
						f.write(str(Name) )
						f.close()
						f1 = open("dict1.txt","w")
						f1.write(str(Rev) )
						f1.close()
					else:
						for i in range(0,y_train.shape[0]):
							y_train.iloc[i] = Name[y_train.iloc[i]] 
					for i in range(0,y_validation.shape[0]):
						y_validation.iloc[i] = Name[y_validation.iloc[i]] 
					### need to convert pandas dataframe and series to np arrays
					X_train = np.asarray(X_train).astype('float32')
					y_train = np.asarray(y_train).astype('int')
					X_validation = np.asarray(X_validation).astype('float32')
					y_validation = np.asarray(y_validation).astype('int')
					model = keras.models.Sequential(
						[keras.layers.Input(shape=X_train.shape[1:]),# input
						keras.layers.Dropout(rate=Dropout_rate),
						keras.layers.Dense(200, activation=activation, 
											kernel_initializer=kernel_initializer,
											kernel_constraint=keras.constraints.max_norm(1.)),
						keras.layers.BatchNormalization(),
						keras.layers.Dropout(rate=Dropout_rate),
						keras.layers.Dense(100, activation=activation, 
											kernel_initializer=kernel_initializer,
											kernel_constraint=keras.constraints.max_norm(1.)),        
						keras.layers.BatchNormalization(),
						keras.layers.Dropout(rate=Dropout_rate),
						keras.layers.Dense(85, activation='softmax')])      # output: 0, 1, 2, 3
					loss = keras.losses.sparse_categorical_crossentropy
					lr_sch = keras.optimizers.schedules.ExponentialDecay(
						initial_learning_rate=1e-2, decay_steps=10000, decay_rate=0.9)
					opti = keras.optimizers.Adam(learning_rate=lr_sch, beta_1=0.9, beta_2=0.999)
					metr = keras.metrics.sparse_categorical_accuracy
					model.compile(loss=loss,      # loss function
								  optimizer=opti, # backpropagation, simple gradient descent 
								  metrics=[metr]) # model performance measure
					checkpoint_path = "./Checkpoint"
					callback_modelcp = keras.callbacks.ModelCheckpoint(filepath=checkpoint_path,save_best_only=True)
					callback_earlystop = keras.callbacks.EarlyStopping(patience=10,restore_best_weights=True)
					callbacks = [callback_modelcp, callback_earlystop]
					history = model.fit(X_train, y_train,epochs=50, verbose=1,validation_data=(X_validation, y_validation),callbacks=callbacks)
					y_pred = model.predict(X_validation)
					y_pred = np.argmax(y_pred, axis=1)
					f1_valid = f1_score(y_validation, y_pred, average='macro')
					if '%s_%s_%s'%(activation, kernel_initializer, Dropout_rate) not in CV_performance:
						CV_performance['%s_%s_%s'%(activation, kernel_initializer, Dropout_rate)] = []
					CV_performance['%s_%s_%s'%(activation, kernel_initializer, Dropout_rate)].append(f1_valid) 

	# get the average F1 for each parameter set
	CV_p = {}
	for parameter in CV_performance:
		if np.mean(CV_performance[parameter]) not in CV_p:
			CV_p[np.mean(CV_performance[parameter])] = []
		CV_p[np.mean(CV_performance[parameter])].append(parameter)

	parameter = CV_p[max(CV_p.keys())][0]
	activation,kernel_initializer,Dropout_rate = parameter.split('_')
	Dropout_rate = float(Dropout_rate)
	
	####### Run ML models #######
	## Make empty dataframes
	conf_matrices = pd.DataFrame(columns = np.insert(arr = classes.astype(np.str), obj = 0, values = 'Class'))
	# imp = pd.DataFrame(index = list(df.drop(['Class'], axis=1)))
	accuracies = []
	accuracies_ho = []
	f1_array = np.array([np.insert(arr = classes.astype(np.str), obj = 0, values = 'M')])
	accuracies_ho = []
	f1_array_ho = np.array([np.insert(arr = test_classes.astype(np.str), obj = 0, values = 'M')])
	results = []
	results_ho = []
	df_all = df.copy()
	df_proba_test = pd.DataFrame(data=df_test['Class'], index=df_test.index, columns=['Class'])
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
		for i in range(0,y_train.shape[0]):
			y_train.iloc[i] = Name[y_train.iloc[i]] 
		for i in range(0,y_validation.shape[0]):
			y_validation.iloc[i] = Name[y_validation.iloc[i]] 
		# for i in range(0,y_test.shape[0]):
			# y_test.iloc[i] = Name[y_test.iloc[i]] 
		### need to convert pandas dataframe and series to np arrays
		X_train = np.asarray(X_train).astype('float32')
		y_train = np.asarray(y_train).astype('int')
		X_validation = np.asarray(X_validation).astype('float32')
		y_validation = np.asarray(y_validation).astype('int')
		X_test = np.asarray(X_test).astype('float32')
		# y_test = np.asarray(y_test).astype('int')
		model = keras.models.Sequential(
			[keras.layers.Input(shape=X_train.shape[1:]),# input
			keras.layers.Dropout(rate=Dropout_rate),
			keras.layers.Dense(200, activation=activation, 
				kernel_initializer=kernel_initializer,
				kernel_constraint=keras.constraints.max_norm(1.)),
			keras.layers.BatchNormalization(),
			keras.layers.Dropout(rate=Dropout_rate),
			keras.layers.Dense(100, activation=activation, 
				kernel_initializer=kernel_initializer,
				kernel_constraint=keras.constraints.max_norm(1.)),        
			keras.layers.BatchNormalization(),
			keras.layers.Dropout(rate=Dropout_rate),
			keras.layers.Dense(85, activation='softmax')])      # output: 0, 1, 2, 3
		loss = keras.losses.sparse_categorical_crossentropy
		lr_sch = keras.optimizers.schedules.ExponentialDecay(
		initial_learning_rate=1e-2, decay_steps=10000, decay_rate=0.9)
		opti = keras.optimizers.Adam(learning_rate=lr_sch, beta_1=0.9, beta_2=0.999)
		metr = keras.metrics.sparse_categorical_accuracy
		model.compile(loss=loss,optimizer=opti,metrics=[metr])
		checkpoint_path = "./Checkpoint"
		callback_modelcp = keras.callbacks.ModelCheckpoint(filepath=checkpoint_path,save_best_only=True)
		callback_earlystop = keras.callbacks.EarlyStopping(patience=10,restore_best_weights=True)
		callbacks = [callback_modelcp, callback_earlystop]
		history = model.fit(X_train, y_train,epochs=50, verbose=1,validation_data=(X_validation, y_validation),callbacks=callbacks)
		y_pred = model.predict(X_validation)
		y_pred = np.argmax(y_pred, axis=1)
		f1_valid = f1_score(y_validation, y_pred, average='macro')
		proba = model.predict(X_validation)
		pred = model.predict_classes(X_validation)
		model.save(save_path + short_name + "_DNN_SMOTE_%s_%s.h5"%(cv_number,dataset))
		pred_va = pd.DataFrame(y_validation)
		pred_va['Prediction'] = pred
		if cv_number == 1:
			Prediction_validation = pred_va
		else:
			Prediction_validation = pd.concat([Prediction_validation,pred_va],axis=0)
		ho_proba = model.predict(X_test)
		ho_pred = model.predict_classes(X_test)
		Prediction_testing['Prediction_%s'%cv_number] = ho_pred
		# convert back to pd series
		y_validation = pd.Series(y_validation)
		pred = pd.Series(pred)
		ho_pred = pd.Series(ho_pred)
		for i in range(0,len(y_validation)):
			y_validation[i] = Rev[y_validation[i]] 
		for i in range(0,len(pred)):
			pred[i] = Rev[pred[i]] 
		for i in range(0,len(ho_pred)):
			ho_pred[i] = Rev[ho_pred[i]] 
		result = Performance_MC(y_validation, pred, classes)
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
		scores_file = save_path + short_name + "_DNN_%s_validation_scores.txt"%dataset
		out_scores = open(scores_file,"w")
		out_scores.write("#ID\t"+pd.DataFrame.to_csv(df_proba,sep="\t").strip()+"\n")
		out_scores.close()
		
		scores_file = save_path + short_name + "_DNN_%s_test_scores.txt"%dataset
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
		out = open(save_path + short_name + "_DNN_%s_results.txt"%dataset, 'w')
		out.write('ID: %s\nAlgorithm: NN\nTrained on classes: %s\nNumber of features: %i\n' % (DF, classes, n_features))
		out.write('Parameters used:activation:%s\nkernel_initializer:%s\nDropout_rate:%s\n' % (activation,kernel_initializer,Dropout_rate))

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
		
		Prediction_validation.to_csv(save_path + short_name + "_DNN_%s_validation_prediction.txt"%dataset,index=True, header=True,sep="\t")
		Prediction_testing.to_csv(save_path + short_name + "_DNN_%s_testing_prediction.txt"%dataset,index=True, header=True,sep="\t")
		

if __name__ == '__main__':
	main()
