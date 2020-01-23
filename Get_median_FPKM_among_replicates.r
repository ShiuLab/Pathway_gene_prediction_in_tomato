setwd('D:\\Pathway_prediction\\Co_expression_20180123')
dat = read.table("Cufflink_matrix_sample_name_20180123.txt", sep="\t", row.names=1,head=T,stringsAsFactors=F)
dat <- dat[,order(colnames(dat))]
#############################################################################
### mutation
### 2012_TFL
res_mutation_FPKM <- c()
subdat <- dat[,21:24]
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}

### 2014_rin_FUL1_2
subdat <- dat[,323:340]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}

### 2015_biosynthesis
subdat <- dat[,351:362]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2015_MSH1
subdat <- dat[,571:582]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2016_SlRBZ
subdat <- dat[,777:784]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2013_P450
subdat <- dat[,113:184]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2014_GLK
subdat <- dat[,198:265]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2015_AGO1
subdat <- dat[,341:344]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2015_ARF3
subdat <- dat[,347:350]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2016_CsMYBF1
subdat <- dat[,700:703]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2016_sulfurea
subdat <- dat[,785:792]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2017_ltm
subdat <- dat[,795:806]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2017_sldml2
subdat <- dat[,869:871]
colnames(subdat)
colnames(subdat)[1:2] <- substr(colnames(subdat)[1:2],1,nchar(colnames(subdat)[1:2])-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2017_ys
subdat <- dat[,896:901]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}


### 2014_leaf_dev
subdat <- dat[,c(266:271,290:292)]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,out)
		colnames(res_mutation_FPKM) <- col_name} else {
		col_name <- c(colnames(res_mutation_FPKM),sam[i])
		res_mutation_FPKM <- cbind(res_mutation_FPKM,tem)
		colnames(res_mutation_FPKM) <- col_name}}

rownames(res_mutation_FPKM) <- rownames(dat)
write.table(res_mutation_FPKM,'Results_median_FPKM_for_mutation_Sly_20180125.txt',quote=F,sep='\t')


#############################################################################
### hormone
### 2013_ABA
res_hormone_FPKM <- c()
subdat <- dat[,25:31]
colnames(subdat)
colnames(subdat)[1] <- substr(colnames(subdat)[1],1,nchar(colnames(subdat)[1])-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,out)
		colnames(res_hormone_FPKM) <- col_name} else {
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,tem)
		colnames(res_hormone_FPKM) <- col_name}}


### 2013_cytokinin
subdat <- dat[,59:76]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-4)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,out)
		colnames(res_hormone_FPKM) <- col_name} else {
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,tem)
		colnames(res_hormone_FPKM) <- col_name}}


### 2013_cytokinin_auxin
subdat <- dat[,77:112]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-4)
colnames(subdat)[c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36)] <- substr(colnames(subdat)[c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36)],1,nchar(colnames(subdat)[c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36)])-1)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,out)
		colnames(res_hormone_FPKM) <- col_name} else {
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,tem)
		colnames(res_hormone_FPKM) <- col_name}}


### 2015_proGRAS
subdat <- dat[,593:600]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,out)
		colnames(res_hormone_FPKM) <- col_name} else {
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,tem)
		colnames(res_hormone_FPKM) <- col_name}}


### 2017_TIBA
subdat <- dat[,872:877]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,out)
		colnames(res_hormone_FPKM) <- col_name} else {
		col_name <- c(colnames(res_hormone_FPKM),sam[i])
		res_hormone_FPKM <- cbind(res_hormone_FPKM,tem)
		colnames(res_hormone_FPKM) <- col_name}}

rownames(res_hormone_FPKM) <- rownames(dat)
write.table(res_hormone_FPKM,'Results_median_FPKM_for_hormone_Sly_20180125.txt',quote=F,sep='\t')

#############################################################################
### environmental, stress, non-stress
### 2013_bacterial
res_stress_FPKM <- c()
subdat <- dat[,32:58]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2014_pathogen
subdat <- dat[,296:322]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2015_cold
subdat <- dat[,448:450]
colnames(subdat)
#colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2015_DC3000
subdat <- dat[,451:462]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2015_fungal
subdat <- dat[,477:479]
colnames(subdat)
#colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2016_A2AS_heat
subdat <- dat[,672:679]
colnames(subdat)
#colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2016_shade_sun
subdat <- dat[,744:776]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2013_TYLCV
subdat <- dat[,185:188]
colnames(subdat)
#colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2015_light
subdat <- dat[,480:570]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2015_pollen_heat
subdat <- dat[,583:592]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2015_T3
subdat <- dat[,613:618]
colnames(subdat)
#colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2016_AvrPto_pathogen
subdat <- dat[,680:699]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2017_fungi
subdat <- dat[,793:794]
colnames(subdat)
#colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2017_pathosystem
subdat <- dat[,807:862]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2017_PSTVd
subdat <- dat[,863:868]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}


### 2015_circadian
subdat <- dat[,363:447]
colnames(subdat)
colnames(subdat) <- gsub('B', 'A', colnames(subdat))
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,out)
		colnames(res_stress_FPKM) <- col_name} else {
		col_name <- c(colnames(res_stress_FPKM),sam[i])
		res_stress_FPKM <- cbind(res_stress_FPKM,tem)
		colnames(res_stress_FPKM) <- col_name}}

rownames(res_stress_FPKM) <- rownames(dat)
write.table(res_stress_FPKM,'Results_median_FPKM_for_stress_Sly_20180125.txt',quote=F,sep='\t')

#############################################################################
### developmental
### 2015_circadian
### see in stress

### 2014_GLK
res_develop_FPKM <- c()
subdat <- dat[,198:265]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2015_fruit
subdat <- dat[,463:476]
colnames(subdat)
#colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2012_fruit_dev
subdat <- dat[,1:20]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2014_leaf_dev
subdat <- dat[,272:295]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2015_anther
subdat <- dat[,345:346]
colnames(subdat)
#colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2016_early_fruit
subdat <- dat[,704:719]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2014_ARF
subdat <- dat[,189:197]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2015_root
subdat <- dat[,601:612]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2015_time_course
subdat <- dat[,619:665]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2015_trichome
subdat <- dat[,666:671]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2016_inflorescence
subdat <- dat[,732:743]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2016_fruit_ripening
subdat <- dat[,720:731]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}


### 2017_yield
subdat <- dat[,878:895]
colnames(subdat)
colnames(subdat) <- substr(colnames(subdat),1,nchar(colnames(subdat))-2)
colnames(subdat)
sam <- unique(colnames(subdat))
for(i in 1:length(sam)){
	tem <- subdat[,colnames(subdat)==sam[i]]
	out <- c()
	if(length(tem) < 100){
		for(j in 1:nrow(tem)){
			out <- c(out,median(as.numeric(tem[j,])))}
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,out)
		colnames(res_develop_FPKM) <- col_name} else {
		col_name <- c(colnames(res_develop_FPKM),sam[i])
		res_develop_FPKM <- cbind(res_develop_FPKM,tem)
		colnames(res_develop_FPKM) <- col_name}}

rownames(res_develop_FPKM) <- rownames(dat)
write.table(res_develop_FPKM,'Results_median_FPKM_for_development_Sly_20180125.txt',quote=F,sep='\t')

res1 <- cbind(res_mutation_FPKM,res_hormone_FPKM)
res1 <- cbind(res1,res_stress_FPKM)
res1 <- cbind(res1,res_develop_FPKM)
write.table(res1,'Results_median_FPKM_for_all_Sly_20180125.txt',quote=F,sep='\t')
