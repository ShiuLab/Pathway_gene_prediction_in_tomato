# Pathway_gene_prediction_in_tomato
Scripts for our paper: __Wang P, Moore BM__, Ugyun S, __Lehti-Shiu M__, Barry C, __Shiu SH__. 	Optimizing the use of gene expression data to predict plant metabolic pathway memberships. *[bioRxiv](https://doi.org/10.1101/2020.07.15.204222)* 


> Expression datasets please see https://zenodo.org/deposit/4581142#
 - Files startswith "Results_" are Set A data
 - Files with "_CV_#_features.txt", where # indicates which cross-validation is (1-5), are Set B data for unsupervised and supervised approaches
 
> The orignal Set B data for naive approaches, and for the Set B used in unsupervised and supervised approaches, please see XXX
 
> ## data preprocessing

> Fold changes
 - Rscript Get_fold_change_values.r

> Median FPKM among replicates
 - Rscript Get_median_FPKM_among_replicates.r

> TPM
 - Rscript Get_TPM.r

> split of genes to test, training and validation (5-fold cross-validation)
 - python Split_data_five_CV.py

> make features for cross-validation sets, Set B
 - python Get_features_for_training_validating_testing.py path_of_expression_matrix expression_matrix number_of_cv

> Get background F1s
 - python Get_background_F1_for_validation_and_test.py

> ## naive approaches
 - python Naive_prediction_crossvalidation_output_prediction.py -path /mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Finished/ -expression Multiclass_MR_pcc_FC_stress_complete.txt

> ## Unsupervised approaches
Unsupervised approaches including:
	kmean
	affinity
	birch
	meanshift
> Example for Set A
 - python Clustering_crossvalidation.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path /mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_nonoverlapping_genes/ -save_path ./Final_results_kmean_setA/ -clustering_method kmean -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> Example for Set B
 - python Clustering_crossvalidation.py -df_short_name spearman_FPKM_hormone -path ./Features_files/ -save_path ./Final_results_kmean_setB/ -clustering_method kmean -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setB
 

> ## Supervised approaches

> RandomForest
 - python RF_crossvalidation.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path /mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_nonoverlapping_genes/ -save_path ./Final_results_RF_setA/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> RandomForest, balance numbers of genes for each pathway for training set
 - python RF_crossvalidation_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path /mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_nonoverlapping_genes/ -save_path ./Final_results_RF_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> SVC, balance numbers of genes for each pathway for training set
 - python SVC_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path /mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_nonoverlapping_genes/ -save_path ./Final_results_SVC_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> KNeighbors, balance numbers of genes for each pathway for training set
 - python KNeighbors_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path /mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_nonoverlapping_genes/ -save_path ./Final_results_SVC_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> KNeighbors, balance numbers of genes for each pathway for training set
 - python KNeighbors_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path /mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_nonoverlapping_genes/ -save_path ./Final_results_KNN_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> DNN, balance numbers of genes for each pathway for training set
 - python DNN_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path /mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_nonoverlapping_genes/ -save_path ./Final_results_DNN_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> neural_network.MLPClassifier, balance numbers of genes for each pathway for training set
 - python neural_network.MLPClassifier_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path /mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_nonoverlapping_genes/ -save_path ./Final_results_MLP_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

