# Pathway_gene_prediction_in_tomato
Scripts for our paper: __Wang P, Moore BM__, Ugyun S, __Lehti-Shiu M__, Barry C, __Shiu SH__. 	Optimizing the use of gene expression data to predict plant metabolic pathway memberships. *[bioRxiv](https://doi.org/10.1101/2020.07.15.204222)* 


> Expression datasets please see https://zenodo.org/record/4585635#.YEKEjMhvrrM
 - File "Cufflink_matrix_sample_name_20180123.txt" contains the original FPKM for each sample replicate
 - File "HTseq_matrix_sample_name_20180123.txt" contains the original read counts for each sample replicate
 - File "Tomato_TPM.txt" contains the original TPM for each sample replicate
 - Files starting with "Results_Fold_changes_", "Results_median_FPKM_" or "TPM_" contain FC, FPKM or TPM values for 41 datasets used in the manuscript.
 
## Get expression matrix

> Fold changes
 - Rscript Get_fold_change_values.r

> Median FPKM among replicates
 - Rscript Get_median_FPKM_among_replicates.r

> Median TPM among replicates
 - Rscript Get_TPM.r

## Calculate the gene-to-gene co-expression matrix

> Partial correlation
 - Rscript Corpcor.r expression_matrix 
##### Note that, if the Coexpression_multiclass.py is used, you don't have to calculate the PCC, Spearman's and Mutual information for the gene-to-gene coexpression matrix using the following three scripts
> Pearson correlation coefficient
 - python  PCC_pandas_20180918.py -file expression_matrix -path path_to_expression_matrix

> Spearman's rank correlation coefficient
 - python Spearman_pandas_20180918.py -file expression_matrix -path path_to_expression_matrix

> Mutual information
 - python MI_pandas_20180919.py -file expression_matrix -path path_to_expression_matrix -start where_the_subset_starts -stop where_the_subset_stops

## Calculate both the gene-to-gene and gene-to-pathway (median or max) co-expression matrix
 - python Coexpression_multiclass.py -pathway_anotation Sly_pathway_annotation_20190117_with_expression_5_members_nonoverlapping.txt -exp Results_median_FPKM_for_stress_Sly_20180125.txt -method pcc

## Calculate the gene-to-gene and gene-to-pathway (median or max) mutual rank of expression similarity
- python Get_MR.py -pathway_anotation Sly_pathway_annotation_20190117_with_expression_5_members_nonoverlapping.txt -exp Results_median_FPKM_for_2017_TIBA.txt -method pcc

## Data preprocessing before model building
> split of genes to test, training and validation (5-fold cross-validation)
 - python Split_data_five_CV.py

> make features for cross-validation sets for Set B
 - python Get_features_for_training_validating_testing.py path_of_expression_matrix expression_matrix number_of_cv

> Get background F1s
 - python Get_background_F1_for_validation_and_test.py

## Naive approaches
 - python Naive_prediction_crossvalidation_output_prediction.py -path ./ -expression Multiclass_MR_pcc_FC_stress.txt

## Unsupervised approaches
Unsupervised approaches including: kmean (KMeans), affinity (AffinityPropagation), birch (Birch), meanshift (MeanShift)
> Example for Set A
 - python Clustering_crossvalidation.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path ./ -save_path ./Final_results_kmean_setA/ -clustering_method kmean -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> Example for Set B
 - python Clustering_crossvalidation.py -df_short_name spearman_FPKM_hormone -path ./Features_files/ -save_path ./Final_results_kmean_setB/ -clustering_method kmean -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setB
 
## Supervised approaches

> RandomForest
 - python RF_crossvalidation.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path ./ -save_path ./Final_results_RF_setA/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> RandomForest, balance numbers of genes for each pathway for training set
 - python RF_crossvalidation_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path ./ -save_path ./Final_results_RF_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> SVC, balance numbers of genes for each pathway for training set
 - python SVC_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path ./ -save_path ./Final_results_SVC_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> KNeighbors, balance numbers of genes for each pathway for training set
 - python KNeighbors_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path ./ -save_path ./Final_results_SVC_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> KNeighbors, balance numbers of genes for each pathway for training set
 - python KNeighbors_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path ./ -save_path ./Final_results_KNN_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> DNN, balance numbers of genes for each pathway for training set
 - python DNN_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path ./ -save_path ./Final_results_DNN_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

> neural_network.MLPClassifier, balance numbers of genes for each pathway for training set
 - python neural_network.MLPClassifier_SMOTE.py -df_short_name Results_median_FPKM_for_stress_Sly_20180125.txt -path ./ -save_path ./Final_results_MLP_setA_SMOTE/ -test_gene_list Genes_for_testing.txt -train_gene_list Genes_for_training.txt -dataset setA

