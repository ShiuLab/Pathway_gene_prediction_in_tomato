setwd('D:\\Pathway_prediction\\Co_expression_20180123')
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("limma")
library("limma")
library("edgeR")
counts = read.table("HTseq_matrix_sample_name_20180123.txt", sep="\t", row.names=1,head=T,stringsAsFactors=F)
counts <- counts[,order(colnames(counts))]
#############################################################################
### mutation
### 2012_TFL
subdat <- counts[,21:24]
group <- factor(c("TFL","TFL","wt","wt"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
FG_de.com <- exactTest(d, pair=c("TFL", "wt"))
results_TFL <- topTags(FG_de.com,n = length(d$AveLogCPM))
res_mutation <- cbind('gene'=rownames(results_TFL$table[1]),'2012_TFL'=results_TFL$table[1])
colnames(res_mutation) <- c('gene','2012_TFL')

### 2014_rin_FUL1_2
subdat <- counts[,323:340]
colnames(subdat)
group <- factor(c("ful_green","ful_green","ful_green","ful_pink","ful_pink","ful_pink","rin_green","rin_green","rin_green","rin_pink","rin_pink","rin_pink","wt_green","wt_green","wt_green","wt_pink","wt_pink","wt_pink")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#ful_green
FG_de.com <- exactTest(d, pair=c("ful_green", "wt_green"))
results_ful_green <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_ful_green$table[1]),results_ful_green$table[1])
colnames(tem) <- c('gene','2014_rin_FUL1_2_ful_green')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
#ful_pink
FG_de.com <- exactTest(d, pair=c("ful_pink", "wt_pink"))
results_ful_pink <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_ful_pink$table[1]),results_ful_pink$table[1])
colnames(tem) <- c('gene','2014_rin_FUL1_2_ful_pink')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
#rin_green
FG_de.com <- exactTest(d, pair=c("rin_green", "wt_green"))
results_rin_green <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_rin_green$table[1]),results_rin_green$table[1])
colnames(tem) <- c('gene','2014_rin_FUL1_2_rin_green')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
#rin_pink
FG_de.com <- exactTest(d, pair=c("rin_pink", "wt_pink"))
results_rin_pink <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_rin_pink$table[1]),results_rin_pink$table[1])
colnames(tem) <- c('gene','2014_rin_FUL1_2_rin_pink')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
### 20170911 comparasion among treatment
treat <- c("ful_green","ful_pink","rin_green","rin_pink")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2014_rin_FUL1_2',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
	}}

### 2015_biosynthesis
subdat <- counts[,351:362]
colnames(subdat)
group <- factor(c("AtMYB12_1week","AtMYB12_1week","Del_Ros_1week","Del_Ros_1week","Del_Ros_4week","Del_Ros_4week","Indigo_1week","Indigo_1week","wt_1week","wt_1week","wt_4week","wt_4week")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#AtMYB12_1week
FG_de.com <- exactTest(d, pair=c("AtMYB12_1week", "wt_1week"))
results_AtMYB12_1week <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_AtMYB12_1week$table[1]),results_AtMYB12_1week$table[1])
colnames(tem) <- c('gene','2015_biosynthesis_AtMYB12_1week')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
#Del_Ros_1week
FG_de.com <- exactTest(d, pair=c("Del_Ros_1week", "wt_1week"))
results_Del_Ros_1week <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_Del_Ros_1week$table[1]),results_Del_Ros_1week$table[1])
colnames(tem) <- c('gene','2015_biosynthesis_Del_Ros_1week')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
#Del_Ros_4week
FG_de.com <- exactTest(d, pair=c("Del_Ros_4week", "wt_4week"))
results_Del_Ros_4week <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_Del_Ros_4week$table[1]),results_Del_Ros_4week$table[1])
colnames(tem) <- c('gene','2015_biosynthesis_Del_Ros_4week')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
#Indigo_1week
FG_de.com <- exactTest(d, pair=c("Indigo_1week", "wt_1week"))
results_Indigo_1week <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_Indigo_1week$table[1]),results_Indigo_1week$table[1])
colnames(tem) <- c('gene','2015_biosynthesis_Indigo_1week')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
### 20170911 comparasion among treatment
treat <- c("AtMYB12_1week","Del_Ros_1week","Del_Ros_4week","Indigo_1week")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_biosynthesis',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
	}}

### 2015_MSH1
subdat <- counts[,571:582]
colnames(subdat)
group <- factor(c("epiF3","epiF3","epiF3","MSH1_Dwarf_dr","MSH1_Dwarf_dr","MSH1_Dwarf_dr","MSH1_mild_dr","MSH1_mild_dr","MSH1_mild_dr","wt","wt","wt")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# epiF3
FG_de.com <- exactTest(d, pair=c("epiF3", "wt"))
results_epiF3 <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_epiF3$table[1]),results_epiF3$table[1])
colnames(tem) <- c('gene','2015_MSH1_epiF3')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# MSH1_Dwarf_dr
FG_de.com <- exactTest(d, pair=c("MSH1_Dwarf_dr", "wt"))
results_MSH1_Dwarf_dr <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_MSH1_Dwarf_dr$table[1]),results_MSH1_Dwarf_dr$table[1])
colnames(tem) <- c('gene','2015_MSH1_MSH1_Dwarf_dr')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# MSH1_mild_dr
FG_de.com <- exactTest(d, pair=c("MSH1_mild_dr", "wt"))
results_MSH1_mild_dr <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_MSH1_mild_dr$table[1]),results_MSH1_mild_dr$table[1])
colnames(tem) <- c('gene','2015_MSH1_MSH1_mild_dr')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
### 20170911 comparasion among treatment
treat <- c("epiF3","MSH1_Dwarf_dr","MSH1_mild_dr")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_MSH1',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
	}}

### 2016_SlRBZ
subdat <- counts[,777:784]
colnames(subdat)
group <- factor(c("SIRBZ","SIRBZ","SIRBZ","SIRBZ","wt","wt","wt","wt")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# SIRBZ
FG_de.com <- exactTest(d, pair=c("SIRBZ", "wt"))
results_SIRBZ <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_SIRBZ$table[1]),results_SIRBZ$table[1])
colnames(tem) <- c('gene','2016_SlRBZ')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')

### 2013_P450
subdat <- counts[,113:184]
colnames(subdat)
group <- factor(c("RNAi_2_SlKLUH_flowerbud","RNAi_2_SlKLUH_flowerbud","RNAi_2_SlKLUH_flowerbud","RNAi_2G2_7dpa_pericarp","RNAi_2G2_7dpa_pericarp","RNAi_2G2_7dpa_pericarp","RNAi_2G2_7dpa_seed","RNAi_2G2_7dpa_seed","RNAi_2G2_7dpa_seed","RNAi_2Q1_7dpa_pericarp","RNAi_2Q1_7dpa_pericarp","RNAi_2Q1_7dpa_pericarp","RNAi_2Q1_7dpa_seed","RNAi_2Q1_7dpa_seed","RNAi_2Q1_7dpa_seed","RNAi_3_SlKLUH_flowerbud","RNAi_3_SlKLUH_flowerbud","RNAi_3_SlKLUH_flowerbud","wt_fw3.2_10dpa_pericarp","wt_fw3.2_10dpa_pericarp","wt_fw3.2_10dpa_pericarp","wt_fw3.2_10dpa_pericarp","wt_fw3.2_10dpa_seed","wt_fw3.2_10dpa_seed","wt_fw3.2_10dpa_seed","wt_fw3.2_10dpa_seed","wt_fw3.2_5dpa_pericarp","wt_fw3.2_5dpa_pericarp","wt_fw3.2_5dpa_pericarp","wt_fw3.2_5dpa_pericarp","wt_fw3.2_5dpa_seed","wt_fw3.2_5dpa_seed","wt_fw3.2_5dpa_seed","wt_fw3.2_5dpa_seed","wt_fw3.2_7dpa_pericarp","wt_fw3.2_7dpa_pericarp","wt_fw3.2_7dpa_pericarp","wt_fw3.2_7dpa_pericarp","wt_fw3.2_7dpa_seed","wt_fw3.2_7dpa_seed","wt_fw3.2_7dpa_seed","wt_fw3.2_7dpa_seed","wt_fw3.2_NIL_flowerbud","wt_fw3.2_NIL_flowerbud","wt_fw3.2_NIL_flowerbud","ys_fw3.2_10dpa_pericarp","ys_fw3.2_10dpa_pericarp","ys_fw3.2_10dpa_pericarp","ys_fw3.2_10dpa_pericarp","ys_fw3.2_10dpa_seed","ys_fw3.2_10dpa_seed","ys_fw3.2_10dpa_seed","ys_fw3.2_10dpa_seed","ys_fw3.2_5dpa_pericarp","ys_fw3.2_5dpa_pericarp","ys_fw3.2_5dpa_pericarp","ys_fw3.2_5dpa_pericarp","ys_fw3.2_5dpa_seed","ys_fw3.2_5dpa_seed","ys_fw3.2_5dpa_seed","ys_fw3.2_5dpa_seed","ys_fw3.2_7dpa_pericarp","ys_fw3.2_7dpa_pericarp","ys_fw3.2_7dpa_pericarp","ys_fw3.2_7dpa_pericarp","ys_fw3.2_7dpa_seed","ys_fw3.2_7dpa_seed","ys_fw3.2_7dpa_seed","ys_fw3.2_7dpa_seed","ys_fw3.2_NIL_flowerbud","ys_fw3.2_NIL_flowerbud","ys_fw3.2_NIL_flowerbud")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# RNAi_2_SlKLUH_flowerbud
FG_de.com <- exactTest(d, pair=c("RNAi_2_SlKLUH_flowerbud", "wt_fw3.2_NIL_flowerbud"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_RNAi_2_SlKLUH_flowerbud')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# RNAi_2G2_7dpa_pericarp
FG_de.com <- exactTest(d, pair=c("RNAi_2G2_7dpa_pericarp", "wt_fw3.2_7dpa_pericarp"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_RNAi_2G2_7dpa_pericarp')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# RNAi_2G2_7dpa_seed
FG_de.com <- exactTest(d, pair=c("RNAi_2G2_7dpa_seed", "wt_fw3.2_7dpa_seed"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_RNAi_2G2_7dpa_seed')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# RNAi_2Q1_7dpa_pericarp
FG_de.com <- exactTest(d, pair=c("RNAi_2Q1_7dpa_pericarp", "wt_fw3.2_7dpa_pericarp"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_RNAi_2Q1_7dpa_pericarp')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# RNAi_2Q1_7dpa_seed
FG_de.com <- exactTest(d, pair=c("RNAi_2Q1_7dpa_seed", "wt_fw3.2_7dpa_seed"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_RNAi_2Q1_7dpa_seed')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# RNAi_3_SlKLUH_flowerbud
FG_de.com <- exactTest(d, pair=c("RNAi_3_SlKLUH_flowerbud", "wt_fw3.2_NIL_flowerbud"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_RNAi_3_SlKLUH_flowerbud')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# ys_fw3.2_10dpa_pericarp
FG_de.com <- exactTest(d, pair=c("ys_fw3.2_10dpa_pericarp", "wt_fw3.2_10dpa_pericarp"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_ys_fw3.2_10dpa_pericarp')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# ys_fw3.2_10dpa_seed
FG_de.com <- exactTest(d, pair=c("ys_fw3.2_10dpa_seed", "wt_fw3.2_10dpa_seed"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_ys_fw3.2_10dpa_seed')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# ys_fw3.2_5dpa_pericarp
FG_de.com <- exactTest(d, pair=c("ys_fw3.2_5dpa_pericarp", "wt_fw3.2_5dpa_pericarp"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_ys_fw3.2_5dpa_pericarp')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# ys_fw3.2_5dpa_seed
FG_de.com <- exactTest(d, pair=c("ys_fw3.2_5dpa_seed", "wt_fw3.2_5dpa_seed"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_ys_fw3.2_5dpa_seed')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# ys_fw3.2_7dpa_pericarp
FG_de.com <- exactTest(d, pair=c("ys_fw3.2_7dpa_pericarp", "wt_fw3.2_7dpa_pericarp"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_ys_fw3.2_7dpa_pericarp')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# ys_fw3.2_7dpa_seed
FG_de.com <- exactTest(d, pair=c("ys_fw3.2_7dpa_seed", "wt_fw3.2_7dpa_seed"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_ys_fw3.2_7dpa_seed')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# ys_fw3.2_NIL_flowerbud
FG_de.com <- exactTest(d, pair=c("ys_fw3.2_NIL_flowerbud", "wt_fw3.2_NIL_flowerbud"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2013_P450_ys_fw3.2_NIL_flowerbud')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
### comparasion among treatment
treat <- c("RNAi_2_SlKLUH_flowerbud","RNAi_2G2_7dpa_pericarp","RNAi_2G2_7dpa_seed","RNAi_2Q1_7dpa_pericarp","RNAi_2Q1_7dpa_seed","RNAi_3_SlKLUH_flowerbud","ys_fw3.2_10dpa_pericarp","ys_fw3.2_10dpa_seed","ys_fw3.2_5dpa_pericarp","ys_fw3.2_5dpa_seed","ys_fw3.2_7dpa_pericarp","ys_fw3.2_7dpa_seed","ys_fw3.2_NIL_flowerbud")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2013_P450',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
	}}

### 2014_GLK
subdat <- counts[,198:265]
colnames(subdat)
group <- factor(c("Breaker_10_bottom","Breaker_10_bottom","Breaker_10_bottom","Breaker_10_middle","Breaker_10_middle","Breaker_10_middle","Breaker_10_top","Breaker_10_top","Breaker_10_top","Breaker_10dpa_bottom","Breaker_10dpa_bottom","Breaker_10dpa_bottom","Breaker_10dpa_control","Breaker_10dpa_control","Breaker_10dpa_middle","Breaker_10dpa_middle","Breaker_10dpa_middle","Breaker_10dpa_OE_GLK1","Breaker_10dpa_OE_GLK1","Breaker_10dpa_OE_GLK2","Breaker_10dpa_OE_GLK2","Breaker_10dpa_OE_GLK2","Breaker_10dpa_top","Breaker_10dpa_top","Breaker_10dpa_top","Breaker_20dpa_bottom","Breaker_20dpa_bottom","Breaker_20dpa_bottom","Breaker_20dpa_middle","Breaker_20dpa_middle","Breaker_20dpa_middle","Breaker_20dpa_top","Breaker_20dpa_top","Breaker_5_bottom","Breaker_5_bottom","Breaker_5_bottom","Breaker_5_middle","Breaker_5_middle","Breaker_5_middle","Breaker_5_top","Breaker_5_top","Breaker_5_top","Breaker_fruit_bottom","Breaker_fruit_bottom","Breaker_fruit_bottom","Breaker_fruit_middle","Breaker_fruit_middle","Breaker_fruit_middle","Breaker_fruit_top","Breaker_fruit_top","Breaker_fruit_top","Breaker_green_fruit_bottom","Breaker_green_fruit_bottom","Breaker_green_fruit_bottom","Breaker_green_fruit_middle","Breaker_green_fruit_middle","Breaker_green_fruit_middle","Breaker_green_fruit_top","Breaker_green_fruit_top","Breaker_green_fruit_top","Immature_green_control","Immature_green_control","Immature_green_control","Immature_green_OE_GLK1","Immature_green_OE_GLK1","Immature_green_OE_GLK2","Immature_green_OE_GLK2","Immature_green_OE_GLK2")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# Breaker_10dpa_OE_GLK1
FG_de.com <- exactTest(d, pair=c("Breaker_10dpa_OE_GLK1", "Breaker_10dpa_control"))
results_epiF3 <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_epiF3$table[1]),results_epiF3$table[1])
colnames(tem) <- c('gene','2014_GLK_Breaker_10dpa_OE_GLK1')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# Breaker_10dpa_OE_GLK2
FG_de.com <- exactTest(d, pair=c("Breaker_10dpa_OE_GLK2", "Breaker_10dpa_control"))
results_epiF3 <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_epiF3$table[1]),results_epiF3$table[1])
colnames(tem) <- c('gene','2014_GLK_Breaker_10dpa_OE_GLK2')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# Immature_green_OE_GLK1
FG_de.com <- exactTest(d, pair=c("Immature_green_OE_GLK1", "Immature_green_control"))
results_epiF3 <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_epiF3$table[1]),results_epiF3$table[1])
colnames(tem) <- c('gene','2014_GLK_Immature_green_OE_GLK1')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# Immature_green_OE_GLK2
FG_de.com <- exactTest(d, pair=c("Immature_green_OE_GLK2", "Immature_green_control"))
results_epiF3 <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_epiF3$table[1]),results_epiF3$table[1])
colnames(tem) <- c('gene','2014_GLK_Immature_green_OE_GLK2')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
### comparasion among treatment
treat <- c("Breaker_10dpa_OE_GLK1","Breaker_10dpa_OE_GLK2","Immature_green_OE_GLK1","Immature_green_OE_GLK2")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2014_GLK',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
	}}

### 2015_AGO1
subdat <- counts[,341:344]
colnames(subdat)
group <- factor(c("AGO1_leaf","AGO1_leaf","TRV_leaf","TRV_leaf")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# AGO1_leaf
FG_de.com <- exactTest(d, pair=c("AGO1_leaf", "TRV_leaf"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2015_AGO1_AGO1_leaf')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')

### 2015_ARF3
subdat <- counts[,347:350]
colnames(subdat)
group <- factor(c("trichome_APF3","trichome_APF3","trichome_WT","trichome_WT")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# trichome_APF3
FG_de.com <- exactTest(d, pair=c("trichome_APF3", "trichome_WT"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2015_ARF3_trichome_APF3')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')

### 2016_CsMYBF1
subdat <- counts[,700:703]
colnames(subdat)
group <- factor(c("Overexpressed_fruit","Overexpressed_fruit","wt_fruit","wt_fruit")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# Overexpressed_fruit
FG_de.com <- exactTest(d, pair=c("Overexpressed_fruit", "wt_fruit"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2016_CsMYBF1_Overexpressed_fruit')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')

### 2016_sulfurea
subdat <- counts[,785:792]
colnames(subdat)
group <- factor(c("sulf_heter","sulf_heter","sulf_heter","sulf_homo","sulf_homo","sulf_homo","wt","wt")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# sulf_heter
FG_de.com <- exactTest(d, pair=c("sulf_heter", "wt"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2016_sulfurea_sulf_heter')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# sulf_homo
FG_de.com <- exactTest(d, pair=c("sulf_homo", "wt"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2016_sulfurea_sulf_homo')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# sulf_homo Vs sulf_heter
FG_de.com <- exactTest(d, pair=c("sulf_homo", "sulf_heter"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2016_sulfurea_sulf_homo_sulf_heter')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')

### 2017_ltm
subdat <- counts[,795:806]
colnames(subdat)
group <- factor(c("ltm_4Leaves","ltm_4Leaves","ltm_6Leaves","ltm_6Leaves","WT_4Leaves","WT_4Leaves","WT_6Leaves","WT_6Leaves","WT_6Leaves","WT_6Leaves","WT_TM","WT_TM")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# ltm_4Leaves
FG_de.com <- exactTest(d, pair=c("ltm_4Leaves", "WT_4Leaves"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2017_ltm_ltm_4Leaves')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# ltm_6Leaves
FG_de.com <- exactTest(d, pair=c("ltm_6Leaves", "WT_6Leaves"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2017_ltm_ltm_6Leaves')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# ltm_6Leaves Vs ltm_4Leaves
FG_de.com <- exactTest(d, pair=c("ltm_6Leaves", "ltm_4Leaves"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2017_ltm_ltm_6Leaves_ltm_4Leaves')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# WT_TM
FG_de.com <- exactTest(d, pair=c("WT_TM", "WT_4Leaves"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2017_ltm_WT_TM')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')

### 2017_sldml2
subdat <- counts[,869:871]
colnames(subdat)
group <- factor(c("sldml2_1_46dpa","sldml2_1_46dpa","WT_46dpa")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# sldml2_1_46dpa
FG_de.com <- exactTest(d, pair=c("sldml2_1_46dpa", "WT_46dpa"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2017_sldml2_sldml2_1_46dpa')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')

### 2017_ys
subdat <- counts[,896:901]
colnames(subdat)
group <- factor(c("stigma_wt","stigma_wt","stigma_wt","stigma_ys","stigma_ys","stigma_ys")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# stigma_ys
FG_de.com <- exactTest(d, pair=c("stigma_ys", "stigma_wt"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2017_ys_stigma_ys')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')

### 2014_leaf_dev
subdat <- counts[,c(266:271,290:292)]
colnames(subdat)
group <- factor(c("35S_BOPa","35S_BOPa","35S_BOPa","BOP_RNAi","BOP_RNAi","BOP_RNAi","wt","wt","wt")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# 35S_BOPa
FG_de.com <- exactTest(d, pair=c("35S_BOPa", "wt"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2014_leaf_dev_35S_BOPa')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# BOP_RNAi
FG_de.com <- exactTest(d, pair=c("BOP_RNAi", "wt"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2014_leaf_dev_BOP_RNAi')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')
# BOP_RNAi Vs 35S_BOPa
FG_de.com <- exactTest(d, pair=c("BOP_RNAi", "35S_BOPa"))
results_tem <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_tem$table[1]),results_tem$table[1])
colnames(tem) <- c('gene','2014_leaf_dev_BOP_RNAi_35S_BOPa')
res_mutation <- merge(res_mutation,tem,by.x='gene',by.y='gene')


#############################################################################
### hormone
### 2013_ABA
subdat <- counts[,25:31]
group <- factor(c("a1d","a2d","c0d","c1d","c2d","ABA_leaf","wt_leaf"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
# ABA_leaf
FG_de.com <- exactTest(d, pair=c("ABA_leaf", "wt_leaf"),dispersion=0.05320461)
results_ABA <- topTags(FG_de.com,n = length(d$AveLogCPM))
res_hormone <- cbind('gene'=rownames(results_ABA$table[1]),'2013_ABA_ABA_leaf'=results_ABA$table[1])
colnames(res_hormone) <- c('gene','2013_ABA')
#a1d
FG_de.com <- exactTest(d, pair=c("a1d", "c0d"),dispersion=0.05320461)
results_auxin_LR <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_auxin_LR$table[1]),results_auxin_LR$table[1])
colnames(tem) <- c('gene','2013_ABA_a1d')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
#a2d
FG_de.com <- exactTest(d, pair=c("a2d", "c2d"),dispersion=0.05320461)
results_auxin_LR <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_auxin_LR$table[1]),results_auxin_LR$table[1])
colnames(tem) <- c('gene','2013_ABA_a2d')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
#a2d Vs a1d
FG_de.com <- exactTest(d, pair=c("a2d", "a1d"),dispersion=0.05320461)
results_auxin_LR <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_auxin_LR$table[1]),results_auxin_LR$table[1])
colnames(tem) <- c('gene','2013_ABA_a2d_a1d')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')

### 2013_cytokinin
subdat <- counts[,59:76]
colnames(subdat)
group <- factor(c("13d_BL24h","13d_BL24h","13d_BL24h","13d_BL2h","13d_BL2h","13d_BL2h","13d_DL24h","13d_DL24h","13d_DL24h","13d_DL2h","13d_DL2h","13d_DL2h","35d_BL24h","35d_BL2h","35d_BL2h","35d_DL24h","35d_DL2h","35d_DL2h"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- c("13d_BL24h","13d_BL2h","35d_BL24h","35d_BL2h","13d_DL24h","13d_DL2h","35d_DL24h","35d_DL2h")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2013_cytokinin',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
	}}

### 2013_cytokinin_auxin
subdat <- counts[,77:112]
colnames(subdat)
group <- factor(c("auxin_LR","auxin_LR","auxin_RT","auxin_RT","auxin_WR","auxin_WR","auxin_LR","auxin_LR","auxin_RT","auxin_RT","auxin_WR","auxin_WR","cytokinin_LR","cytokinin_LR","cytokinin_RT","cytokinin_RT","cytokinin_WR","cytokinin_WR","cytokinin_LR","cytokinin_LR","cytokinin_RT","cytokinin_RT","cytokinin_WR","cytokinin_WR","DMSO_LR","DMSO_LR","DMSO_LR","DMSO_LR","DMSO_RT","DMSO_RT","DMSO_RT","DMSO_RT","DMSO_WR","DMSO_WR","DMSO_WR","DMSO_WR"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#auxin_LR
FG_de.com <- exactTest(d, pair=c("auxin_LR", "DMSO_LR"))
results_auxin_LR <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_auxin_LR$table[1]),results_auxin_LR$table[1])
colnames(tem) <- c('gene','2013_cytokinin_auxin_auxin_LR')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
#auxin_RT
FG_de.com <- exactTest(d, pair=c("auxin_RT", "DMSO_RT"))
results_auxin_RT <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_auxin_RT$table[1]),results_auxin_RT$table[1])
colnames(tem) <- c('gene','2013_cytokinin_auxin_auxin_RT')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
#auxin_WR
FG_de.com <- exactTest(d, pair=c("auxin_WR", "DMSO_WR"))
results_auxin_WR <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_auxin_WR$table[1]),results_auxin_WR$table[1])
colnames(tem) <- c('gene','2013_cytokinin_auxin_auxin_WR')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
#cytokinin_LR
FG_de.com <- exactTest(d, pair=c("cytokinin_LR", "DMSO_LR"))
results_cytokinin_LR <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_cytokinin_LR$table[1]),results_cytokinin_LR$table[1])
colnames(tem) <- c('gene','2013_cytokinin_auxin_cytokinin_LR')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
#cytokinin_RT
FG_de.com <- exactTest(d, pair=c("cytokinin_RT", "DMSO_RT"))
results_cytokinin_RT <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_cytokinin_RT$table[1]),results_cytokinin_RT$table[1])
colnames(tem) <- c('gene','2013_cytokinin_auxin_cytokinin_RT')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
#cytokinin_WR
FG_de.com <- exactTest(d, pair=c("cytokinin_WR", "DMSO_WR"))
results_cytokinin_WR <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_cytokinin_WR$table[1]),results_cytokinin_WR$table[1])
colnames(tem) <- c('gene','2013_cytokinin_auxin_cytokinin_WR')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
### 20170911 comparasion among treatment
treat <- c("auxin_LR","auxin_RT","auxin_WR","cytokinin_LR","cytokinin_RT","cytokinin_WR")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2013_cytokinin_auxin',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
	}}

### 2015_proGRAS
subdat <- counts[,593:600]
colnames(subdat)
group <- factor(c("proGRAS_Pac_GA","proGRAS_Pac_GA","proGRAS_Pac","proGRAS_Pac","wt_Pac_GA","wt_Pac_GA","wt_Pac","wt_Pac")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# wt_Pac_GA
FG_de.com <- exactTest(d, pair=c("wt_Pac_GA", "wt_Pac"))
results_wt_Pac_GA <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_Pac_GA$table[1]),results_wt_Pac_GA$table[1])
colnames(tem) <- c('gene','2015_proGRAS_wt_Pac_GA')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
# proGRAS_Pac
FG_de.com <- exactTest(d, pair=c("proGRAS_Pac", "wt_Pac"))
results_proGRAS_Pac <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_proGRAS_Pac$table[1]),results_proGRAS_Pac$table[1])
colnames(tem) <- c('gene','2015_proGRAS_proGRAS_Pac')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
treat <- c("wt_Pac_GA","proGRAS_Pac","proGRAS_Pac_GA")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_proGRAS',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
	}}

### 2017_TIBA
subdat <- counts[,872:877]
colnames(subdat)
group <- factor(c("lfs","lfs","TIBA_pin","TIBA_pin","wt","wt")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# lfs
FG_de.com <- exactTest(d, pair=c("lfs", "wt"))
results_wt_Pac_GA <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_Pac_GA$table[1]),results_wt_Pac_GA$table[1])
colnames(tem) <- c('gene','2017_TIBA_lfs')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
# TIBA_pin
FG_de.com <- exactTest(d, pair=c("TIBA_pin", "wt"))
results_wt_Pac_GA <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_Pac_GA$table[1]),results_wt_Pac_GA$table[1])
colnames(tem) <- c('gene','2017_TIBA_pin')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')
# TIBA_pin Vs lfs
FG_de.com <- exactTest(d, pair=c("TIBA_pin", "lfs"))
results_wt_Pac_GA <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_Pac_GA$table[1]),results_wt_Pac_GA$table[1])
colnames(tem) <- c('gene','2017_TIBA_pin_lfs')
res_hormone <- merge(res_hormone,tem,by.x='gene',by.y='gene')


#############################################################################
### environmental, stress, non-stress
### 2013_bacterial
subdat <- counts[,32:58]
colnames(subdat)
group <- factor(c("Agrobacterium_tumefaciens","Agrobacterium_tumefaciens","Agrobacterium_tumefaciens","bacteria_mock","bacteria_mock","bacteria_mock","DC3000","DC3000","DC3000","DC3000AvrPtoAvrPtoB","DC3000AvrPtoAvrPtoB","DC3000AvrPtoAvrPtoB","DC3000hrcQ_UfliC","DC3000hrcQ_UfliC","DC3000hrcQ_UfliC","flgII28","flgII28","flgII28","PAMPs_mock","PAMPs_mock","PAMPs_mock","Pseudomonas_fluorescens","Pseudomonas_fluorescens","Pseudomonas_fluorescens","Pseudomonas_putida","Pseudomonas_putida","Pseudomonas_putida"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#flgII28
FG_de.com <- exactTest(d, pair=c("flgII28", "PAMPs_mock"))
results_flgII28 <- topTags(FG_de.com,n = length(d$AveLogCPM))
res_stress <- cbind(rownames(results_flgII28$table[1]),results_flgII28$table[1])
colnames(res_stress) <- c('gene','2013_bacterial_flgII28')
### 20170911 comparasion with mock
mock <- "bacteria_mock"
treat <- c("Agrobacterium_tumefaciens","DC3000","DC3000AvrPtoAvrPtoB","DC3000hrcQ_UfliC","Pseudomonas_fluorescens","Pseudomonas_putida")
for(i in 1:length(treat)){
	FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), "bacteria_mock"))
	results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
	tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
	colnames(tem) <- c('gene',paste('2013_bacterial_',as.character(treat[i]),sep=''))
	res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}
### 20170911 comparasion among treatment
treat <- c("Agrobacterium_tumefaciens","DC3000","DC3000AvrPtoAvrPtoB","DC3000hrcQ_UfliC","flgII28","Pseudomonas_fluorescens","Pseudomonas_putida")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2013_bacterial',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}


### 2014_pathogen
subdat <- counts[,296:322]
colnames(subdat)
group <- factor(c("prf19_DC3000_4h","prf19_DC3000_4h","prf19_DC3000_4h","prf19_DC3000_6h","prf19_DC3000_6h","prf19_DC3000_6h","prf3_DC3000_4h","prf3_DC3000_4h","prf3_DC3000_4h","prf3_DC3000_6h","prf3_DC3000_6h","prf3_DC3000_6h","PtoR_DC3000_4h","PtoR_DC3000_4h","PtoR_DC3000_4h","PtoR_DC3000_6h","PtoR_DC3000_6h","PtoR_DC3000_6h","PtoR_DC3000_B","PtoR_DC3000_B","PtoR_DC3000_B","PtoR_DC3000_C","PtoR_DC3000_C","PtoR_DC3000_C","PtoR_DC3000_D","PtoR_DC3000_D","PtoR_DC3000_D")) ### PtoR_DC3000_B, PtoR_DC3000_C, PtoR_DC3000_D were not used
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#prf19_DC3000_4h
FG_de.com <- exactTest(d, pair=c("prf19_DC3000_4h", "PtoR_DC3000_4h"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2014_pathogen_prf19_DC3000_4h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#prf19_DC3000_6h
FG_de.com <- exactTest(d, pair=c("prf19_DC3000_6h", "PtoR_DC3000_6h"))
results_prf19_DC3000_6h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_6h$table[1]),results_prf19_DC3000_6h$table[1])
colnames(tem) <- c('gene','2014_pathogen_prf19_DC3000_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#prf3_DC3000_4h
FG_de.com <- exactTest(d, pair=c("prf3_DC3000_4h", "PtoR_DC3000_4h"))
results_prf3_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf3_DC3000_4h$table[1]),results_prf3_DC3000_4h$table[1])
colnames(tem) <- c('gene','2014_pathogen_prf3_DC3000_4h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#prf3_DC3000_6h
FG_de.com <- exactTest(d, pair=c("prf3_DC3000_6h", "PtoR_DC3000_6h"))
results_prf3_DC3000_6h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf3_DC3000_6h$table[1]),results_prf3_DC3000_6h$table[1])
colnames(tem) <- c('gene','2014_pathogen_prf3_DC3000_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
### 20170911 comparasion among treatment
treat <- c("prf19_DC3000_4h","prf19_DC3000_6h","prf3_DC3000_4h","prf3_DC3000_6h")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2014_pathogen',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}

### 2015_cold
subdat <- counts[,448:450]
colnames(subdat)
group <- factor(c("0h","12h","1h")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
# 1h
FG_de.com <- exactTest(d, pair=c("1h", "0h"),dispersion=0.05320461)
results_1h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind('gene'=rownames(results_1h$table[1]),results_1h$table[1])
colnames(tem) <- c('gene','2015_cold_1h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# 12h
FG_de.com <- exactTest(d, pair=c("12h", "0h"),dispersion=0.05320461)
results_12h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind('gene'=rownames(results_12h$table[1]),results_12h$table[1])
colnames(tem) <- c('gene','2015_cold_12h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# 1h Vs 12h
FG_de.com <- exactTest(d, pair=c("1h", "12h"),dispersion=0.05320461)
results_12h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind('gene'=rownames(results_12h$table[1]),results_12h$table[1])
colnames(tem) <- c('gene','2015_cold_1h_12h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')


### 2015_DC3000
subdat <- counts[,451:462]
colnames(subdat)
group <- factor(c("DC3000_leaf","DC3000_leaf","DC3000_leaf","leaf","leaf","leaf","red_light_DC3000_leaf","red_light_DC3000_leaf","red_light_DC3000_leaf","red_light_leaf","red_light_leaf","red_light_leaf")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- c("DC3000_leaf","red_light_DC3000_leaf","red_light_leaf", "leaf")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_DC3000',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}

### 2015_fungal
subdat <- counts[,477:479]
colnames(subdat)
group <- factor(c("Appressoria","Quiescent","wt")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
### 20170911 comparasion among treatment
treat <- c("Appressoria","Quiescent", "wt")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])),dispersion=0.05320461)
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_fungal',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}

### 2016_A2AS_heat
subdat <- counts[,672:679]
colnames(subdat)
group <- factor(c("A2AS_anther_control","A2AS_anther_heat","A2AS_leaf_control","A2AS_leaf_heat","wt_anther_control","wt_anther_heat","wt_leaf_control","wt_leaf_heat")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
# A2AS_anther
FG_de.com <- exactTest(d, pair=c("A2AS_anther_control", "wt_anther_control"),dispersion=0.05320461)
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2016_A2AS_heat_A2AS_anther')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# A2AS_anther_heat
FG_de.com <- exactTest(d, pair=c("A2AS_anther_heat", "A2AS_anther_control"),dispersion=0.05320461)
results_A2AS_anther_heat <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther_heat$table[1]),results_A2AS_anther_heat$table[1])
colnames(tem) <- c('gene','2016_A2AS_heat_A2AS_anther_heat')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# A2AS_leaf
FG_de.com <- exactTest(d, pair=c("A2AS_leaf_control", "wt_leaf_control"),dispersion=0.05320461)
results_A2AS_leaf <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_leaf$table[1]),results_A2AS_leaf$table[1])
colnames(tem) <- c('gene','2016_A2AS_heat_A2AS_leaf')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# A2AS_leaf_heat
FG_de.com <- exactTest(d, pair=c("A2AS_leaf_heat", "A2AS_leaf_control"),dispersion=0.05320461)
results_A2AS_leaf_heat <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_leaf_heat$table[1]),results_A2AS_leaf_heat$table[1])
colnames(tem) <- c('gene','2016_A2AS_heat_A2AS_leaf_heat')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# wt_anther_heat
FG_de.com <- exactTest(d, pair=c("wt_anther_heat", "wt_anther_control"),dispersion=0.05320461)
results_wt_anther_heat <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_anther_heat$table[1]),results_wt_anther_heat$table[1])
colnames(tem) <- c('gene','2016_A2AS_heat_wt_anther_heat')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# wt_leaf_heat
FG_de.com <- exactTest(d, pair=c("wt_leaf_heat", "wt_leaf_control"),dispersion=0.05320461)
results_wt_leaf_heat <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_leaf_heat$table[1]),results_wt_leaf_heat$table[1])
colnames(tem) <- c('gene','2016_A2AS_heat_wt_leaf_heat')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
### 20170911 comparasion among treatment
treat <- c("A2AS_anther_heat","A2AS_leaf_heat","wt_anther_heat","wt_leaf_heat")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])),dispersion=0.05320461)
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2016_A2AS_heat',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}

### 2016_shade_sun
subdat <- counts[,744:776]
colnames(subdat)
group <- factor(c("seedling","seedling","seedling","seedling","seedling","seedling","seedling","shade_floral","shade_floral","shade_floral","shade_leaf","shade_seedling","shade_seedling","shade_shoot","shade_shoot","shade_stem","shade_stem","shade_veg_meristrem","shade_veg_meristrem","shade_veg_meristrem","sun_floral","sun_floral","sun_leaf","sun_leaf","sun_seedling","sun_seedling","sun_shoot","sun_shoot","sun_shoot","sun_stem","sun_stem","sun_veg_meristrem","sun_veg_meristrem")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
# shade_seedling
FG_de.com <- exactTest(d, pair=c("shade_seedling", "seedling"))
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2016_shade_sun_shade_seedling')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# sun_seedling
FG_de.com <- exactTest(d, pair=c("sun_seedling", "seedling"))
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2016_shade_sun_sun_seedling')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# shade_seedling
FG_de.com <- exactTest(d, pair=c("shade_seedling", "sun_seedling"))
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2016_shade_sun_seedling')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# shade_floral
FG_de.com <- exactTest(d, pair=c("shade_floral", "sun_floral"))
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2016_shade_sun_floral')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# shade_leaf
FG_de.com <- exactTest(d, pair=c("shade_leaf", "sun_leaf"))
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2016_shade_sun_leaf')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# shade_shoot
FG_de.com <- exactTest(d, pair=c("shade_shoot", "sun_shoot"))
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2016_shade_sun_shoot')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# shade_stem
FG_de.com <- exactTest(d, pair=c("shade_stem", "sun_stem"))
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2016_shade_sun_stem')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# shade_veg_meristrem
FG_de.com <- exactTest(d, pair=c("shade_veg_meristrem", "sun_veg_meristrem"))
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2016_shade_sun_veg_meristrem')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')

### 2013_TYLCV
subdat <- counts[,185:188]
colnames(subdat)
group <- factor(c("CLN2777A_leaves_control","CLN2777A_leaves_TYLCV","TMXA48_4_0_leaves_control","TMXA48_4_0_leaves_TYLCV")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
# CLN2777A_leaves_TYLCV
FG_de.com <- exactTest(d, pair=c("CLN2777A_leaves_TYLCV", "CLN2777A_leaves_control"),dispersion=0.05320461)
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2013_TYLCV_CLN2777A_leaves_TYLCV')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# TMXA48_4_0_leaves_TYLCV
FG_de.com <- exactTest(d, pair=c("TMXA48_4_0_leaves_TYLCV", "TMXA48_4_0_leaves_control"),dispersion=0.05320461)
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2013_TYLCV_TMXA48_4_0_leaves_TYLCV')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# TMXA48_4_0_leaves_TYLCV Vs CLN2777A_leaves_TYLCV 
FG_de.com <- exactTest(d, pair=c("TMXA48_4_0_leaves_TYLCV", "CLN2777A_leaves_TYLCV"),dispersion=0.05320461)
results_A2AS_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_A2AS_anther$table[1]),results_A2AS_anther$table[1])
colnames(tem) <- c('gene','2013_TYLCV_TMXA48_4_0_leaves_TYLCV_CLN2777A_leaves_TYLCV')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')

### 2015_light
subdat <- counts[,480:570]
colnames(subdat)
group <- factor(c("leaf_primordium_11d","leaf_primordium_11d","leaf_primordium_11d","leaf_primordium_11d","leaf_primordium_11d","leaf_primordium_11d","leaf_primordium_11d","leaf_primordium_14d","leaf_primordium_14d","leaf_primordium_14d","leaf_primordium_14d","leaf_primordium_17d","leaf_primordium_17d","leaf_primordium_17d","leaf_primordium_4d","leaf_primordium_4d","leaf_primordium_4d","leaf_primordium_5d","leaf_primordium_5d","leaf_primordium_5d","leaf_primordium_5d","leaf_primordium_constant_shade","leaf_primordium_constant_shade","leaf_primordium_constant_shade","leaf_primordium_constant_shade","leaf_primordium_constant_shade","leaf_primordium_constant_sun","leaf_primordium_constant_sun","leaf_primordium_constant_sun","leaf_primordium_constant_sun","leaf_primordium_transient_shift","leaf_primordium_transient_shift","leaf_primordium_transient_shift","leaf_primordium_transient_shift","leaf_primordium_transient_shift","leaf_primordium_transient_shift","leaf_primordium_transient_shift","leaf_primordium_transient_shift_shade","leaf_primordium_transient_shift_shade","leaf_primordium_transient_shift_shade","leaf_primordium_transient_shift_shade","leaf_primordium_transient_shift_shade","leaf_primordium_transient_shift_shade","leaf_primordium_transient_shift_shade","SAM_11d","SAM_11d","SAM_11d","SAM_11d","SAM_11d","SAM_11d","SAM_11d","SAM_11d","SAM_14d","SAM_14d","SAM_14d","SAM_14d","SAM_17d","SAM_17d","SAM_17d","SAM_17d","SAM_4d","SAM_4d","SAM_4d","SAM_4d","SAM_5d","SAM_5d","SAM_5d","SAM_5d","SAM_constant_shade","SAM_constant_shade","SAM_constant_shade","SAM_constant_shade","SAM_constant_shade","SAM_constant_sun","SAM_constant_sun","SAM_constant_sun","SAM_constant_sun","SAM_constant_sun","SAM_transient_shift","SAM_transient_shift","SAM_transient_shift","SAM_transient_shift","SAM_transient_shift","SAM_transient_shift","SAM_transient_shift","SAM_transient_shift_shade","SAM_transient_shift_shade","SAM_transient_shift_shade","SAM_transient_shift_shade","SAM_transient_shift_shade","SAM_transient_shift_shade")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20180124 comparasion among treatment
treat <- c("leaf_primordium_11d","leaf_primordium_14d","leaf_primordium_17d","leaf_primordium_4d","leaf_primordium_5d","leaf_primordium_constant_shade","leaf_primordium_constant_sun","leaf_primordium_transient_shift","leaf_primordium_transient_shift_shade")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_light',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}
### 20180124 comparasion among treatment
treat <- c("SAM_11d","SAM_14d","SAM_17d","SAM_4d","SAM_5d","SAM_constant_shade","SAM_constant_sun","SAM_transient_shift","SAM_transient_shift_shade")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_light',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}

### 2015_pollen_heat
subdat <- counts[,583:592]
colnames(subdat)
group <- factor(c("Pol_Control","Pol_Control","Pol_Control","Pol_Control","Pol_Control","Pol_Heat","Pol_Heat","Pol_Heat","Pol_Heat","Pol_Heat")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#Pol_Heat
FG_de.com <- exactTest(d, pair=c("Pol_Heat", "Pol_Control"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2015_pollen_heat')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')

### 2015_T3
subdat <- counts[,613:618]
colnames(subdat)
group <- factor(c("OH88119_mock","OH88119_T3_6d","OH88119_T3_6h","PI114490_mock","PI114490_T3_6d","PI114490_T3_6h")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
# OH88119_T3_6d
FG_de.com <- exactTest(d, pair=c("OH88119_T3_6d", "OH88119_mock"),dispersion=0.05320461)
results_wt_leaf_heat <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_leaf_heat$table[1]),results_wt_leaf_heat$table[1])
colnames(tem) <- c('gene','2015_T3_OH88119_T3_6d')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# OH88119_T3_6h
FG_de.com <- exactTest(d, pair=c("OH88119_T3_6h", "OH88119_mock"),dispersion=0.05320461)
results_wt_leaf_heat <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_leaf_heat$table[1]),results_wt_leaf_heat$table[1])
colnames(tem) <- c('gene','2015_T3_OH88119_T3_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# PI114490_T3_6d
FG_de.com <- exactTest(d, pair=c("PI114490_T3_6d", "PI114490_mock"),dispersion=0.05320461)
results_wt_leaf_heat <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_leaf_heat$table[1]),results_wt_leaf_heat$table[1])
colnames(tem) <- c('gene','2015_T3_PI114490_T3_6d')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
# PI114490_T3_6h
FG_de.com <- exactTest(d, pair=c("PI114490_T3_6h", "PI114490_mock"),dispersion=0.05320461)
results_wt_leaf_heat <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_leaf_heat$table[1]),results_wt_leaf_heat$table[1])
colnames(tem) <- c('gene','2015_T3_PI114490_T3_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
### 20180124 comparasion among treatment
treat <- c("OH88119_T3_6d","OH88119_T3_6h","PI114490_T3_6d","PI114490_T3_6h")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])),dispersion=0.05320461)
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_T3',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}

### 2016_AvrPto_pathogen
subdat <- counts[,680:699]
colnames(subdat)
group <- factor(c("avrPto","avrPto","avrPto_2xA","avrPto_2xA","avrPto_2xA","avrPto_2xA","avrPto","avrPto","avrPto_I96A","avrPto_I96A","avrPto_I96A_2xA","avrPto_I96A_2xA","avrPto_I96A_2xA","avrPto_I96A_2xA","avrPto_I96A","avrPto_I96A","mock","mock","mock","mock")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#avrPto
FG_de.com <- exactTest(d, pair=c("avrPto", "mock"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2016_AvrPto_pathogen_avrPto')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#avrPto_2xA
FG_de.com <- exactTest(d, pair=c("avrPto_2xA", "mock"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2016_AvrPto_pathogen_avrPto_2xA')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#avrPto_I96A
FG_de.com <- exactTest(d, pair=c("avrPto_I96A", "mock"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2016_AvrPto_pathogen_avrPto_I96A')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#avrPto_I96A_2xA
FG_de.com <- exactTest(d, pair=c("avrPto_I96A_2xA", "mock"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2016_AvrPto_pathogen_avrPto_I96A_2xA')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
### 20180124 comparasion among treatment
treat <- c("avrPto","avrPto_2xA","avrPto_I96A","avrPto_I96A_2xA")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2016_AvrPto_pathogen',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}

### 2017_fungi
subdat <- counts[,793:794]
colnames(subdat)
group <- factor(c("control","Infected")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
# Infected
FG_de.com <- exactTest(d, pair=c("Infected", "control"),dispersion=0.05320461)
results_wt_leaf_heat <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_wt_leaf_heat$table[1]),results_wt_leaf_heat$table[1])
colnames(tem) <- c('gene','2017_fungi_Infected')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')

### 2017_pathosystem
subdat <- counts[,807:862]
colnames(subdat)
group <- factor(c("DC300028_6h","DC300028_6h","DC300028_6h","DC3000avrPto_DC3000avrPtoB_4h","DC3000avrPto_DC3000avrPtoB_4h","DC3000avrPto_DC3000avrPtoB_4h","DC3000avrPto_DC3000avrPtoB_6h","DC3000avrPto_DC3000avrPtoB_6h","DC3000avrPto_DC3000avrPtoB_6h","DC3000hopQ1_1_6h","DC3000hopQ1_1_6h","DC3000hopQ1_1_6h","DC3000hrcQ_U_6h","DC3000hrcQ_U_6h","DC3000hrcQ_U_6h","control_bacterial_30m","control_bacterial_30m","control_bacterial_30m","control_bacterial_6h","control_bacterial_6h","control_bacterial_6h","control_PAMPs_30m","control_PAMPs_30m","control_PAMPs_30m","control_PAMPs_6h","control_PAMPs_6h","control_PAMPs_6h","csp22_30m","csp22_30m","csp22_30m","csp22_6h","csp22_6h","csp22_6h","flg22_30m","flg22_30m","flg22_30m","flg22_6h","flg22_6h","flgII_28_30m","flgII_28_30m","flgII_28_30m","K12_30m","K12_30m","K12_30m","K12_6h","K12_6h","K12_6h","MgCl2_30m","MgCl2_30m","MgCl2_30m","SigmaL2262_30m","SigmaL2262_30m","SigmaL2262_30m","SigmaL2262_6h","SigmaL2262_6h","SigmaL2262_6h")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#DC300028_6h
FG_de.com <- exactTest(d, pair=c("DC300028_6h", "control_bacterial_6h"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_DC300028_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#DC3000avrPto_DC3000avrPtoB_4h
FG_de.com <- exactTest(d, pair=c("DC3000avrPto_DC3000avrPtoB_4h", "control_bacterial_6h"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_DC3000avrPto_DC3000avrPtoB_4h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#DC3000avrPto_DC3000avrPtoB_6h
FG_de.com <- exactTest(d, pair=c("DC3000avrPto_DC3000avrPtoB_6h", "control_bacterial_6h"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_DC3000avrPto_DC3000avrPtoB_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#DC3000hopQ1_1_6h
FG_de.com <- exactTest(d, pair=c("DC3000hopQ1_1_6h", "control_bacterial_6h"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_DC3000hopQ1_1_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#DC3000hrcQ_U_6h
FG_de.com <- exactTest(d, pair=c("DC3000hrcQ_U_6h", "control_bacterial_6h"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_DC3000hrcQ_U_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#csp22_30m
FG_de.com <- exactTest(d, pair=c("csp22_30m", "MgCl2_30m"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_csp22_30m')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#csp22_6h
FG_de.com <- exactTest(d, pair=c("csp22_6h", "MgCl2_30m"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_csp22_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#flg22_30m
FG_de.com <- exactTest(d, pair=c("flg22_30m", "control_PAMPs_30m"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_flg22_30m')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#flg22_6h
FG_de.com <- exactTest(d, pair=c("flg22_6h", "control_PAMPs_6h"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_flg22_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#flgII_28_30m
FG_de.com <- exactTest(d, pair=c("flgII_28_30m", "control_PAMPs_30m"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_flgII_28_30m')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#K12_30m
FG_de.com <- exactTest(d, pair=c("K12_30m", "control_PAMPs_30m"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_K12_30m')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#K12_6h
FG_de.com <- exactTest(d, pair=c("K12_6h", "control_PAMPs_6h"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_K12_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#SigmaL2262_30m
FG_de.com <- exactTest(d, pair=c("SigmaL2262_30m", "control_PAMPs_30m"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_SigmaL2262_30m')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#SigmaL2262_6h
FG_de.com <- exactTest(d, pair=c("SigmaL2262_6h", "control_PAMPs_6h"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_pathosystem_SigmaL2262_6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
### 20180124 comparasion among treatment
treat <- c("DC300028_6h","DC3000avrPto_DC3000avrPtoB_4h","DC3000avrPto_DC3000avrPtoB_6h","DC3000hopQ1_1_6h","DC3000hrcQ_U_6h","csp22_30m","csp22_6h","flg22_30m","flg22_6h","flgII_28_30m","K12_30m","K12_6h","SigmaL2262_30m","SigmaL2262_6h")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2017_pathosystem',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}

### 2017_PSTVd
subdat <- counts[,863:868]
colnames(subdat)
group <- factor(c("mock","mock","mock","PSTVd","PSTVd","PSTVd")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#PSTVd
FG_de.com <- exactTest(d, pair=c("PSTVd", "mock"))
results_prf19_DC3000_4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_prf19_DC3000_4h$table[1]),results_prf19_DC3000_4h$table[1])
colnames(tem) <- c('gene','2017_PSTVd')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')

### 2015_circadian
subdat <- counts[,363:447]
colnames(subdat)
group <- factor(c("0h","10h","12h","14h","16h","18h","2h","20h","22h","24h","26h","28h","30h","32h","34h","36h","38h","4h","40h","42h","44h","46h","48h","6h","8h","0h","10h","12h","14h","16h","18h","2h","20h","22h","24h","26h","28h","30h","32h","34h","36h","38h","4h","40h","42h","44h","46h","48h","6h","8h","LL0h","LL10h","LL12h","LL14h","LL16h","LL2h","LL20h","LL22h","LL24h","LL26h","LL28h","LL30h","LL32h","LL34h","LL36h","LL38h","LL4h","LL40h","LL42h","LL44h","LL46h","LL48h","LL50h","LL52h","LL54h","LL56h","LL58h","LL6h","LL60h","LL62h","LL64h","LL66h","LL68h","LL70h","LL8h")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
#LL0h
FG_de.com <- exactTest(d, pair=c("LL0h", "0h"))
results_LL0h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL0h$table[1]),results_LL0h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL0h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL2h
FG_de.com <- exactTest(d, pair=c("LL2h", "2h"))
results_LL2h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL2h$table[1]),results_LL2h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL2h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL4h
FG_de.com <- exactTest(d, pair=c("LL4h", "4h"))
results_LL4h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL4h$table[1]),results_LL4h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL4h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL6h
FG_de.com <- exactTest(d, pair=c("LL6h", "6h"))
results_LL6h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL6h$table[1]),results_LL6h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL6h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL8h
FG_de.com <- exactTest(d, pair=c("LL8h", "8h"))
results_LL8h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL8h$table[1]),results_LL8h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL8h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL10h
FG_de.com <- exactTest(d, pair=c("LL10h", "10h"))
results_LL10h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL10h$table[1]),results_LL10h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL10h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL12h
FG_de.com <- exactTest(d, pair=c("LL12h", "12h"))
results_LL12h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL12h$table[1]),results_LL12h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL12h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL14h
FG_de.com <- exactTest(d, pair=c("LL14h", "14h"))
results_LL14h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL14h$table[1]),results_LL14h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL14h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL16h
FG_de.com <- exactTest(d, pair=c("LL16h", "16h"))
results_LL16h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL16h$table[1]),results_LL16h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL16h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL20h
FG_de.com <- exactTest(d, pair=c("LL20h", "20h"))
results_LL20h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL20h$table[1]),results_LL20h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL20h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL22h
FG_de.com <- exactTest(d, pair=c("LL22h", "22h"))
results_LL22h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL22h$table[1]),results_LL22h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL22h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL24h
FG_de.com <- exactTest(d, pair=c("LL24h", "24h"))
results_LL24h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL24h$table[1]),results_LL24h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL24h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL26h
FG_de.com <- exactTest(d, pair=c("LL26h", "26h"))
results_LL26h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL26h$table[1]),results_LL26h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL26h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL28h
FG_de.com <- exactTest(d, pair=c("LL28h", "28h"))
results_LL28h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL28h$table[1]),results_LL28h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL28h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL30h
FG_de.com <- exactTest(d, pair=c("LL30h", "30h"))
results_LL30h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL30h$table[1]),results_LL30h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL30h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL32h
FG_de.com <- exactTest(d, pair=c("LL32h", "32h"))
results_LL32h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL32h$table[1]),results_LL32h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL32h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL34h
FG_de.com <- exactTest(d, pair=c("LL34h", "34h"))
results_LL34h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL34h$table[1]),results_LL34h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL34h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL36h
FG_de.com <- exactTest(d, pair=c("LL36h", "36h"))
results_LL36h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL36h$table[1]),results_LL36h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL36h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL38h
FG_de.com <- exactTest(d, pair=c("LL38h", "38h"))
results_LL38h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL38h$table[1]),results_LL38h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL38h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL40h
FG_de.com <- exactTest(d, pair=c("LL40h", "40h"))
results_LL40h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL40h$table[1]),results_LL40h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL40h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL42h
FG_de.com <- exactTest(d, pair=c("LL42h", "42h"))
results_LL42h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL42h$table[1]),results_LL42h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL42h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL44h
FG_de.com <- exactTest(d, pair=c("LL44h", "44h"))
results_LL44h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL44h$table[1]),results_LL44h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL44h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL46h
FG_de.com <- exactTest(d, pair=c("LL46h", "46h"))
results_LL46h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL46h$table[1]),results_LL46h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL46h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
#LL48h
FG_de.com <- exactTest(d, pair=c("LL48h", "46h"))
results_LL48h <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_LL48h$table[1]),results_LL48h$table[1])
colnames(tem) <- c('gene','2015_circadian_LL48h')
res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
### 20170911 comparasion among treatment
treat <- c("LL0h","LL2h","LL4h","LL6h","LL8h","LL10h","LL12h","LL14h","LL16h","LL20h","LL22h","LL24h","LL26h","LL28h","LL30h","LL32h","LL34h","LL36h","LL38h","LL40h","LL42h","LL44h","LL46h","LL48h")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_circadian',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_stress <- merge(res_stress,tem,by.x='gene',by.y='gene')
	}}


#############################################################################
### developmental
### 2015_circadian
subdat <- counts[,363:447]
colnames(subdat)
group <- factor(c("0h","10h","12h","14h","16h","18h","2h","20h","22h","24h","26h","28h","30h","32h","34h","36h","38h","4h","40h","42h","44h","46h","48h","6h","8h","0h","10h","12h","14h","16h","18h","2h","20h","22h","24h","26h","28h","30h","32h","34h","36h","38h","4h","40h","42h","44h","46h","48h","6h","8h","LL0h","LL10h","LL12h","LL14h","LL16h","LL2h","LL20h","LL22h","LL24h","LL26h","LL28h","LL30h","LL32h","LL34h","LL36h","LL38h","LL4h","LL40h","LL42h","LL44h","LL46h","LL48h","LL50h","LL52h","LL54h","LL56h","LL58h","LL6h","LL60h","LL62h","LL64h","LL66h","LL68h","LL70h","LL8h")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
res_develop <- c()
### 20180124 comparasion among treatment
treat <- c("0h","10h","12h","14h","16h","18h","2h","20h","22h","24h","26h","28h","30h","32h","34h","36h","38h","4h","40h","42h","44h","46h","48h","6h","8h")
nn = 0
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_circadian',as.character(treat[i]),as.character(treat[j]),sep='_'))
		if(nn==0) {
			res_develop <- tem
			nn = nn + 1}
		else {res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')}
	}}

### 2014_GLK
subdat <- counts[,198:265]
colnames(subdat)
group <- factor(c("Breaker_10_bottom","Breaker_10_bottom","Breaker_10_bottom","Breaker_10_middle","Breaker_10_middle","Breaker_10_middle","Breaker_10_top","Breaker_10_top","Breaker_10_top","Breaker_10dpa_bottom","Breaker_10dpa_bottom","Breaker_10dpa_bottom","Breaker_10dpa_control","Breaker_10dpa_control","Breaker_10dpa_middle","Breaker_10dpa_middle","Breaker_10dpa_middle","Breaker_10dpa_OE_GLK1","Breaker_10dpa_OE_GLK1","Breaker_10dpa_OE_GLK2","Breaker_10dpa_OE_GLK2","Breaker_10dpa_OE_GLK2","Breaker_10dpa_top","Breaker_10dpa_top","Breaker_10dpa_top","Breaker_20dpa_bottom","Breaker_20dpa_bottom","Breaker_20dpa_bottom","Breaker_20dpa_middle","Breaker_20dpa_middle","Breaker_20dpa_middle","Breaker_20dpa_top","Breaker_20dpa_top","Breaker_5_bottom","Breaker_5_bottom","Breaker_5_bottom","Breaker_5_middle","Breaker_5_middle","Breaker_5_middle","Breaker_5_top","Breaker_5_top","Breaker_5_top","Breaker_fruit_bottom","Breaker_fruit_bottom","Breaker_fruit_bottom","Breaker_fruit_middle","Breaker_fruit_middle","Breaker_fruit_middle","Breaker_fruit_top","Breaker_fruit_top","Breaker_fruit_top","Breaker_green_fruit_bottom","Breaker_green_fruit_bottom","Breaker_green_fruit_bottom","Breaker_green_fruit_middle","Breaker_green_fruit_middle","Breaker_green_fruit_middle","Breaker_green_fruit_top","Breaker_green_fruit_top","Breaker_green_fruit_top","Immature_green_control","Immature_green_control","Immature_green_control","Immature_green_OE_GLK1","Immature_green_OE_GLK1","Immature_green_OE_GLK2","Immature_green_OE_GLK2","Immature_green_OE_GLK2")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
treat <- c("Breaker_10_bottom","Breaker_10_middle","Breaker_10_top","Breaker_10dpa_bottom","Breaker_10dpa_middle","Breaker_10dpa_top","Breaker_20dpa_bottom","Breaker_20dpa_middle","Breaker_20dpa_top","Breaker_5_bottom","Breaker_5_middle","Breaker_5_top","Breaker_fruit_bottom","Breaker_fruit_middle","Breaker_fruit_top","Breaker_green_fruit_bottom","Breaker_green_fruit_middle","Breaker_green_fruit_top")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2014_GLK',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2015_fruit
subdat <- counts[,463:476]
colnames(subdat)
group <- factor(c("ac_14daf","ac_21daf","ac_28daf","ac_35daf","ac_42daf","ac_49daf","ac_7daf","HG6_61_14daf","HG6_61_21daf","HG6_61_28daf","HG6_61_35daf","HG6_61_42daf","HG6_61_49daf","HG6_61_7daf")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
### 20170911 comparasion among treatment
treat <- c("ac_7daf","ac_14daf","ac_21daf","ac_28daf","ac_35daf","ac_42daf","ac_49daf","HG6_61_7daf","HG6_61_14daf","HG6_61_21daf","HG6_61_28daf","HG6_61_35daf","HG6_61_42daf","HG6_61_49daf")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])),dispersion=0.05320461)
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_fruit',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2012_fruit_dev
subdat <- counts[,1:20]
group <- factor(c("10d_fruit","10d_fruit","1cm_fruit","1cm_fruit","2cm_fruit","2cm_fruit","3cm_fruit","3cm_fruit","breaker_fruit","breaker_fruit","flower","flower","flower_bud","flower_bud","leaf","leaf","mature_green_fruit","mature_green_fruit","root","root"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- c("10d_fruit","1cm_fruit","2cm_fruit","3cm_fruit","breaker_fruit","flower","flower_bud","leaf","mature_green_fruit","root")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2012_fruit_dev',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2014_leaf_dev
subdat <- counts[,272:295]
group <- factor(c("distal_P4","distal_P4","distal_P4","distal_P5","distal_P5","distal_P5","distal_P6","distal_P6","distal_P6","proximal_P4","proximal_P4","proximal_P4","proximal_P5","proximal_P5","proximal_P5","proximal_P6","proximal_P6","proximal_P6","proximal_P7","proximal_P7","proximal_P7","SAM_P3","SAM_P3","SAM_P3"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- c("distal_P4","distal_P5","distal_P6","proximal_P4","proximal_P5","proximal_P6","proximal_P7","SAM_P3")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2014_leaf_dev',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2015_anther
subdat <- counts[,345:346]
group <- factor(c("anther","leaf")) 
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
# anther
FG_de.com <- exactTest(d, pair=c("anther", "leaf"),dispersion=0.05320461)
results_anther <- topTags(FG_de.com,n = length(d$AveLogCPM))
tem <- cbind(rownames(results_anther$table[1]),results_anther$table[1])
colnames(tem) <- c('gene','2015_anther_anther')
res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')

### 2016_early_fruit
subdat <- counts[,704:719]
group <- factor(c("ov_0dpa","ov_0dpa","ov_1dpa","ov_1dpa","ov_2dpa","ov_2dpa","ov_5dpa","ov_5dpa","ow_0dpa","ow_0dpa","ow_1dpa","ow_1dpa","ow_2dpa","ow_2dpa","ow_5dpa","ow_5dpa"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- c("ov_0dpa","ov_1dpa","ov_2dpa","ov_5dpa","ow_0dpa","ow_1dpa","ow_2dpa","ow_5dpa")
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2016_early_fruit',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2014_ARF
subdat <- counts[,189:197]
group <- factor(c("flower_athesis","flower_athesis","flower_athesis","flower_bud","flower_bud","flower_bud","fruit_4d","fruit_4d","fruit_4d"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- unique(group)
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2014_ARF',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2015_root
subdat <- counts[,601:612]
group <- factor(c("root_hair","root_hair","root_hair","SlDZ","SlDZ","SlDZ","SlEZ","SlEZ","SlEZ","SlM","SlM","SlM"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- unique(group)
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_root',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2015_time_course
subdat <- counts[,619:665]
group <- factor(c("tp12_dark","tp12_dark","tp16_dark","tp16_dark","tp20_dark","tp20_dark","tp24_light","tp24_light","tp28_light","tp28_light","tp32_light","tp32_light","tp36_light","tp36_light","tp40_light","tp40_light","tp44_light","tp44_light","tp48_light","tp48_light","tp52_light","tp52_light","tp56_light","tp56_light","tp60_light","tp60_light","tp64_light","tp64_light","tp64_light","tp68_light","tp68_light","tp68_light","tp72_light","tp72_light","tp72_light","tp72_light","tp76_light","tp76_light","tp76_light","tp76_light","tp80_light","tp80_light","tp80_light","tp80_light","tp84_light","tp84_light","tp84_light"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- unique(group)
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_time_course',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2015_trichome
subdat <- counts[,666:671]
group <- factor(c("shaved_stem","shaved_stem","stem_trichome","stem_trichome","trichome","trichome"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- unique(group)
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2015_trichome',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2016_inflorescence
subdat <- counts[,732:743]
group <- factor(c("EVM","EVM","FM","FM","LVM","LVM","MVM","MVM","SIM","SIM","TM","TM"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- unique(group)
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2016_inflorescence',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2016_fruit_ripening
subdat <- counts[,720:731]
group <- factor(c("fruit_20dpa","fruit_20dpa","fruit_20dpa","fruit_40dpa","fruit_40dpa","fruit_40dpa","fruit_break","fruit_break","fruit_break","fruit_green","fruit_green","fruit_green"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- unique(group)
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2016_fruit_ripening',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}

### 2017_yield
subdat <- counts[,878:895]
group <- factor(c("s_mutant_FM","s_mutant_FM","s_mutant_LVM","s_mutant_LVM","s_mutant_MVM","s_mutant_MVM","s_mutant_SIM","s_mutant_SIM","s_mutant_TM","s_mutant_TM","s2_mutant_FM","s2_mutant_FM","s2_mutant_MVM","s2_mutant_MVM","s2_mutant_SIM","s2_mutant_SIM","s2_mutant_TM","s2_mutant_TM"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d$common.dispersion
### 20170911 comparasion among treatment
treat <- unique(group)
for(i in 1:(length(treat)-1)){
	for(j in (i+1):length(treat)){
		FG_de.com <- exactTest(d, pair=c(as.character(treat[i]), as.character(treat[j])))
		results_treat <- topTags(FG_de.com,n = length(d$AveLogCPM))
		tem <- cbind(rownames(results_treat$table[1]),results_treat$table[1])
		colnames(tem) <- c('gene',paste('2017_yield',as.character(treat[i]),as.character(treat[j]),sep='_'))
		res_develop <- merge(res_develop,tem,by.x='gene',by.y='gene')
	}}


### finally, write the result
write.table(res_mutation,'Results_Fold_changes_for_mutation_Sly_20180124.txt',quote=F,sep='\t',row.names=F)
write.table(res_hormone,'Results_Fold_changes_for_hormone_Sly_20180124.txt',quote=F,sep='\t',row.names=F)
write.table(res_stress,'Results_Fold_changes_for_stress_Sly_20180124.txt',quote=F,sep='\t',row.names=F)
write.table(res_develop,'Results_Fold_changes_for_developmental_Sly_20180124.txt',quote=F,sep='\t',row.names=F)
### combine all the FC together
res_total <- merge(res_develop,res_mutation,by.x='gene',by.y='gene')
res_total_2 <- merge(res_total,res_stress,by.x='gene',by.y='gene')
res_total_3 <- merge(res_total_2,res_hormone,by.x='gene',by.y='gene')
write.table(res_total_3,'Results_Fold_changes_all_combination_Sly_20180124.txt',quote=F,sep='\t',row.names=F)

