ceRNA.Lin <-
function(miRtar,targetce=NULL,geneexp,miRexp,numMIR=1,method="pearson",numrandom=100){
miRtar_d<-split(miRtar[,1],miRtar[,2])
miR_l<-sapply(miRtar_d,length)
miRs_name<-names(miRtar_d)
#ceRNA_O_1<-combinations(length(miRs_name),2,miRs_name,repeats.allowed=F)
t(combn(miRs_name,2))-> ceRNA_O_1
ceRNA<-function(miRtar,targetce=NULL){
if(is.null(targetce)){
#ceRNA_o<-permutations(length(miRs_name),2,miRs_name,repeats.allowed=T)

fun<-function(x){
tmp<-intersect(miRtar_d[[x[1]]],miRtar_d[[x[2]]])
tmp_p<-paste(tmp,collapse=",")
tmp_l<-length(tmp)
return(c(tmp_p,tmp_l))
}#for fun
miRs<-t(apply(ceRNA_O_1,1,fun))
ceRNA<-data.frame(ceRNA_O_1,miRs[,1],as.numeric(miRs[,2]),stringsAsFactors=FALSE)
rownames(ceRNA)<-NULL
colnames(ceRNA)<-c("targetce","anotherce","miRNAs","miRNAs_num")
}else{
tarmiRs<-miRtar_d[[targetce]]
fun1<-function(x){
tmp<-intersect(miRtar_d[[x]],tarmiRs)
tmp_p<-paste(tmp,collapse=",")
tmp_l<-length(tmp)
return(c(tmp_p,tmp_l))
}#for fun1
miRs<-t(sapply(miRs_name,fun1))
ceRNA<-data.frame(targetce,rownames(miRs),miRs[,1],as.numeric(miRs[,2]),stringsAsFactors=FALSE)
rownames(ceRNA)<-NULL
colnames(ceRNA)<-c("targetce","anotherce","miRNAs","miRNAs_num")
}#for else

result<-list(ceRNA=ceRNA,miR_l=miR_l)
return(result)
}
geneexp<-as.matrix(geneexp)
miRexp<-as.matrix(miRexp)
miRtar_geneexp <- intersect(rownames(geneexp),miRtar[,2])
miRtar_miRexp <- intersect(rownames(miRexp),miRtar[,1])
index1 <- miRtar[,2]%in%miRtar_geneexp
miRtar_new <- miRtar[index1,]
index2 <- miRtar_new[,1]%in%miRtar_miRexp
miRtar <- miRtar_new[index2,]

miRtar_d<-split(miRtar[,1],miRtar[,2])
miR_l<-sapply(miRtar_d,length)
#S random 
random.s <- function(mirnum){
s_random <- matrix(nrow=0,ncol=mirnum)
for (i in 1:numrandom){
gene_sample <- geneexp[sample(1:nrow(geneexp),2,replace=FALSE),]
miRexp_index<-sample(1:nrow(miRexp),mirnum,replace=FALSE)
miRNA_sample <- miRexp[miRexp_index,]
gene_names <- rownames(gene_sample)
mir_names <- rownames(miRexp)[miRexp_index]
gene_names <- t(as.matrix(gene_names))
mir_names <- as.matrix(mir_names)
ce <- cbind(matrix(rep(gene_names,each=nrow(mir_names)),nrow=nrow(mir_names)),mir_names)
s <- apply(ce,1,function(j){fun.prepartial(x=geneexp[j[1],],y=geneexp[j[2],],z=miRexp[j[3],])}[2])
s_random <- rbind(s_random,s)
}
s_mean <- apply(s_random,1,function(x){mean(x)})
s_mean <- as.matrix(s_mean)
return(s_mean)
}
#Partial correlation formula
fun.prepartial<-function(x,y,z){
tmp1 <- cor(x,y)
tmp2 <- cor(x,z)
tmp3 <- cor(z,y)
tmp4 <- tmp1-tmp2*tmp3
tmp5 <- sqrt(1-tmp2*tmp2)*sqrt(1-tmp3*tmp3)
partial_cor <- tmp4/tmp5
s_partial <- tmp1-partial_cor
return(c(partial_cor,s_partial))
}#for  fun.prepartial
#Defining fun.partial to calculate the pvalue of Sensitivity value
fun.partial<-function(x){
tmp<-intersect(miRtar_d[[x[1]]],miRtar_d[[x[2]]]) 
tmp_p<-paste(tmp,collapse=",")
tmp_l<-length(tmp)

s_real <- matrix(nrow=0,ncol=1)
for(i in 1:tmp_l){
s_tmp<-fun.prepartial(x=geneexp[x[1],],y=geneexp[x[2],],z=miRexp[tmp[i],])[2]
s_real <- rbind(s_real,s_tmp)
}
sreal_mean <- mean(s_real)
miRnum_background <- random_background[,which(colnames(random_background)==tmp_l)]
pvalue <- length(which(miRnum_background>=sreal_mean))/numrandom

return(c(tmp_p,tmp_l,pvalue))
}#for fun.partial
fun.pearson<-function(x){
tmp<-intersect(miRtar_d[[x[1]]],miRtar_d[[x[2]]])
tmp_p<-paste(tmp,collapse=",")
tmp_l<-length(tmp)
pvalue<-cor.test(geneexp[x[1],],geneexp[x[2],])$p.value
return(c(tmp_p,tmp_l,pvalue))
}#for fun.pearson
if(is.null(targetce)){
ceall <- ceRNA(miRtar=miRtar)$ceRNA
ceall<-ceall[(as.numeric(ceall[,4]))>=numMIR,]
#ceall <- ceall[-which(ceall[,4]==0),]
ceRNA_all <- ceall[,c(1,2)]
miRnum <- unique(ceall[,4])
miRnum <- as.matrix(miRnum)
random_background <- sapply(miRnum,random.s)
colnames(random_background) <- miRnum


if(method=="pearson"){
miRs<-t(apply(ceRNA_all,1,fun.pearson))
}else if(method=="partial correlation"){
miRs<-t(apply(ceRNA_all,1,fun.partial))
}

ceRNA<-data.frame(ceRNA_all,miRs[,1],as.numeric(miRs[,2]),as.numeric(miRs[,3]),stringsAsFactors=FALSE)
#ceRNA_1 <- as.matrix(do.call(rbind,lapply(lapply(seq_len(dim(ceRNA_O_1)[1]), function(i) ceRNA[intersect(which(ceRNA[,1]==ceRNA_O_1[i,1]),which(ceRNA[,2]==ceRNA_O_1[i,2])),]),data.frame)))
rownames(ceRNA)<-NULL
colnames(ceRNA)<-c("targetce","anotherce","miRNAs","miRNAs_num","pvalue")
}else{ 

ce_target <- ceRNA(miRtar=miRtar,targetce=targetce)$ceRNA
ce_target<-ce_target[(as.numeric(ce_target[,4]))>=numMIR,]
#ce_target <- ce_target[-which(ce_target[,4]==0),]
cetarget <- ce_target[,c(1,2)]

miRnum <- unique(ce_target[,4])
miRnum <- as.matrix(miRnum)
random_background <- sapply(miRnum,random.s)
colnames(random_background) <- miRnum

if(method=="pearson"){
miRs<-t(apply(cetarget,1,fun.pearson))
}else if(method=="partial correlation"){


miRs<-t(apply(cetarget,1,fun.partial))
}

ceRNA<-data.frame(cetarget,miRs[,1],as.numeric(miRs[,2]),as.numeric(miRs[,3]),stringsAsFactors=FALSE)
#ceRNA_1 <- as.matrix(do.call(rbind,lapply(lapply(seq_len(dim(ceRNA_O_1)[1]), function(i) ceRNA[intersect(which(ceRNA[,1]==ceRNA_O_1[i,1]),which(ceRNA[,2]==ceRNA_O_1[i,2])),]),data.frame)))
rownames(ceRNA)<-NULL
colnames(ceRNA)<-c("targetce","anotherce","miRNAs","miRNAs_num","pvalue")
}#for else

result<-list(ceRNA=ceRNA,miR_l=miR_l)
return(result)
}
