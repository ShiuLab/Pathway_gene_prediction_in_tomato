setwd('D:\\Pathway_prediction\\20180831\\Within_between_pathways')
df <- read.table('All_within_between_pcc_FPKM_all.txt_02',head=T,sep='\t',stringsAsFactors=F)
df <- df[df[,3]!="PWY-5048",]
new <- df[df[,4]=='within',]
pathway = unique(new[,3])
res_max <- c() ### max PCC for each gene in pathways
res_median <- c() ### median PCC for each gene in pathways
for(i in 1:length(pathway)){
	subdat <- new[new[,3]==pathway[i],]
	gene = unique(c(subdat[,1],subdat[,2]))
	for(j in 1:length(gene)){
		subgene = rbind(subdat[subdat[,1]==gene[j],],subdat[subdat[,2]==gene[j],])
		res_max <- rbind(res_max,c(gene[j],pathway[i],max(subgene[,5])))
		res_median <- rbind(res_median,c(gene[j],pathway[i],median(subgene[,5])))
	}}
colnames(res_max) <- c('Gene','Pathway','Max_PCC')
colnames(res_median) <- c('Gene','Pathway','Max_PCC')
write.table(res_max,'Max_PCC_for_each_gene_within_pathway_20190409.txt',row.names=F,sep='\t',quote=F)
write.table(res_median,'Median_PCC_for_each_gene_within_pathway_20190409.txt',row.names=F,sep='\t',quote=F)
dat_max <- read.table('Median_PCC_for_each_gene_within_pathway_20190409.txt',head=T,sep='\t',stringsAsFactors=F)
dat_median <- read.table('Max_PCC_for_each_gene_within_pathway_20190409.txt',head=T,sep='\t',stringsAsFactors=F)
pathway <- unique(dat_max[,2])
res_max_median <- c()
for(i in 1:length(pathway)){
	res_max_median <- rbind(res_max_median,c(pathway[i],median(dat_max[dat_max[,2]==pathway[i],3])))
	}
res_max_median <- res_max_median[order(as.numeric(res_max_median[,2])),]

res_median_median <- c()
for(i in 1:length(pathway)){
	res_median_median <- rbind(res_median_median,c(pathway[i],median(dat_median[dat_median[,2]==pathway[i],3])))
	}
res_median_median <- res_median_median[order(as.numeric(res_median_median[,2])),]

pcc <- df[df[,4]=='between',] ### 4888351 in total
#pcc <- pcc[sample(nrow(pcc),size=100000,replace=FALSE),]
library(vioplot)
pdf('Vioplot_random_10000_between_pathway_PCC.pdf')
vioplot(pcc[,5])
dev.off()
gene <- unique(c(as.character(pcc[,1]),as.character(pcc[,2])))
pcc_max <- c() ### max PCC for genes in random gene sets
pcc_median <- c() ### max PCC for genes in random gene sets
for(i in 1:length(gene)){
	subdat <- pcc[pcc[,1]==gene[i] | pcc[,2]==gene[i],]
	p <- unique(subdat[,3])
	for(j in 1:length(p)){
		tem <- subdat[subdat[,3]==p[j],]
			# pcc_max <- rbind(pcc_max,c(gene[i],p[j],max(tem[,5])))
			# pcc_median <- rbind(pcc_median,c(gene[i],p[j],median(tem[,5])))
			pcc_max <- c(pcc_max,max(tem[,5]))
			pcc_median <- c(pcc_median,median(tem[,5]))
		}
	}
pdf('Vioplot_random_10000_between_pathway_PCC_max_and_median_for_gene_to_pathways.pdf')
par(mfrow=c(1,2))
vioplot(pcc_median)
vioplot(pcc_max)
dev.off()

ecdf_fun <- function(value_list,value) ecdf(value_list)(value)
per_PCC <- merge(res_median_median,res_max_median,by.x='V1',by.y='V1')
# per_PCC <- cbind(per_PCC,ecdf_fun(as.numeric(as.character(pcc_median[,3])),as.numeric(as.character(per_PCC[,2]))))
# per_PCC <- cbind(per_PCC,ecdf_fun(as.numeric(as.character(pcc_max[,3])),as.numeric(as.character(per_PCC[,3]))))
per_PCC <- cbind(per_PCC,ecdf_fun(as.numeric(as.character(pcc_median)),as.numeric(as.character(per_PCC[,2]))))
per_PCC <- cbind(per_PCC,ecdf_fun(as.numeric(as.character(pcc_max)),as.numeric(as.character(per_PCC[,3]))))
pdf('Median_Max_PCC_and_percentile_within_pathway_20200123.pdf')
plot(as.numeric(as.character(per_PCC[,4])),as.numeric(as.character(per_PCC[,5])),xlim=c(0,1),ylim=c(0,1))
abline(0,1)
dev.off()
write.table(per_PCC,'Median_Max_PCC_and_percentile_within_pathway_20200123.txt',row.names=F,col.names=F,quote=F,sep='\t')	
