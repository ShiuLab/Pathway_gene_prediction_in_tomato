#BiocManager::install("scater")
library(scater)
library(stringr)
library(devtools)
library(Biobase)
library(preprocessCore)
setwd('D:\\Pathway_prediction\\Revision_New_Phytologist\\TMP_calling')
Expression <- read.table('HTseq_matrix_sample_name_20180123_with_length.txt',head=T,sep='\t',stringsAsFactors=F,row.names=1)
Len <- Expression$Length
tpm_table <- calculateTPM(as.matrix(Expression[,1:(ncol(Expression)-1)]), Len)
write.table(tpm_table,'Tomato_TPM.txt',row.names=T,sep='\t',quote=F)

# For the chilling datasets
### tissue_list <- unique(str_remove_all(as.character(colnames(Expression)), "\\.\\d$"))
# median TPM
for(i in seq(0,48,2)){
	colnames(tpm_table) <- gsub(paste('2015_circadian.wt_leaf_B',i,sep=''),paste('2015_circadian.wt_leaf_',i,sep=''),colnames(tpm_table))
	colnames(tpm_table) <- gsub(paste('2015_circadian.wt_leaf_A',i,sep=''),paste('2015_circadian.wt_leaf_',i,sep=''),colnames(tpm_table))
	}
colnames(tpm_table)[which(colnames(tpm_table)=='X2015_circadian.wt_leaf_0')] <- 'X2015_circadian.wt_leaf_0h'
colnames(tpm_table)[which(colnames(tpm_table)=='X2015_circadian.wt_leaf_2')] <- 'X2015_circadian.wt_leaf_2h'
colnames(tpm_table)[which(colnames(tpm_table)=='X2015_circadian.wt_leaf_4')] <- 'X2015_circadian.wt_leaf_4h'
colnames(tpm_table)[which(colnames(tpm_table)=='X2015_circadian.wt_leaf_6')] <- 'X2015_circadian.wt_leaf_6h'
colnames(tpm_table)[which(colnames(tpm_table)=='X2015_circadian.wt_leaf_8')] <- 'X2015_circadian.wt_leaf_8h'
colnames(tpm_table) <- gsub('DL2h_A','DL2h',colnames(tpm_table))
colnames(tpm_table) <- gsub('DL2h_B','DL2h',colnames(tpm_table))
colnames(tpm_table) <- gsub('DL24h_A','DL24h',colnames(tpm_table))
colnames(tpm_table) <- gsub('DL24h_B','DL24h',colnames(tpm_table))
colnames(tpm_table) <- gsub('BL2h_A','BL2h',colnames(tpm_table))
colnames(tpm_table) <- gsub('BL2h_B','BL2h',colnames(tpm_table))
colnames(tpm_table) <- gsub('BL24h_A','BL24h',colnames(tpm_table))
colnames(tpm_table) <- gsub('BL24h_B','BL24h',colnames(tpm_table))
colnames(tpm_table) <- gsub('_WR_II','_WR',colnames(tpm_table))
colnames(tpm_table) <- gsub('_WR_I','_WR',colnames(tpm_table))
colnames(tpm_table) <- gsub('_LR_II','_LR',colnames(tpm_table))
colnames(tpm_table) <- gsub('_LR_I','_LR',colnames(tpm_table))
colnames(tpm_table) <- gsub('_RT_II','_RT',colnames(tpm_table))
colnames(tpm_table) <- gsub('_RT_I','_RT',colnames(tpm_table))
tissue_list <- unique(str_remove_all(as.character(colnames(tpm_table)), "_\\d$"))
tbl_median_counts <- data.frame(row.names = rownames(tpm_table))
for (i in 1:length(tissue_list)) {
  subdat <- tpm_table[,grep(colnames(tpm_table),pattern=tissue_list[i],fixed = TRUE)]
  if(!is.null(nrow(subdat))){
	  new_dat <- as.data.frame(rowMedians(subdat), optional = T)
	  colnames(new_dat) = paste(tissue_list[i],sep='')
	  tbl_median_counts <- cbind(tbl_median_counts, new_dat)
	}
else{
	tbl_median_counts <- cbind(tbl_median_counts, subdat)
	colnames(tbl_median_counts)[i] = paste(tissue_list[i],sep='')
	}
}
tbl_median_counts <- tbl_median_counts[,order(colnames(tbl_median_counts))]
write.table(tbl_median_counts,'TPM_all.txt',sep='\t',quote=F)

colnames(tbl_median_counts) <- gsub('X2013_cytokinin\\.','X2013_cytokinin_.',colnames(tbl_median_counts))
pos <- regexpr('\\.',colnames(tbl_median_counts))
study <- unique(substr(colnames(tbl_median_counts),1,pos-1))
for(i in 1:length(study)){
	subdat <- tbl_median_counts[,grep(colnames(tbl_median_counts),pattern=study[i],fixed = TRUE)]
	if(ncol(subdat)>=3){
		write.table(subdat,paste('TPM_',study[i],'.txt',sep=''),sep='\t',quote=F)
		}
	}

