ceRNA.basic <-
function(miRtar,targetce=NULL,method="ratio",numMIR=1,cutoff=ifelse(method=="ratio",1/3,0.05)){
miRtar_d<-split(miRtar[,1],miRtar[,2])
miR_l<-sapply(miRtar_d,length)

miRs_name<-names(miRtar_d)
#ceRNA_O_1<-combinations(length(miRs_name),2,miRs_name,repeats.allowed=F)
t(combn(miRs_name,2))-> ceRNA_O_1
ceRNA<-function(miRtar,targetce=NULL){
if(is.null(targetce)){
#ceRNA_o<-permutations(length(miRs_name),2,miRs_name,repeats.allowed=T)######permutations

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
if(is.null(targetce)){
tmp<-ceRNA(miRtar=miRtar)
}else{
tmp<-ceRNA(miRtar=miRtar,targetce=targetce)
}#for else
tmpceRNA<-tmp$ceRNA[(as.numeric(tmp$ceRNA[,4]))>=numMIR,]
tmpl_shared<-tmpceRNA[,4]
tmpl_miRs1<-tmp$miR_l[tmpceRNA[,1]]
tmpl_miRs2<-tmp$miR_l[tmpceRNA[,2]]

if(method=="ratio"){
ratio<-tmpl_shared/tmpl_miRs2
ind<-ratio>cutoff
ind1<-!ind
cesig<-data.frame(tmpceRNA[ind,],ratio[ind])
cenotsig<-data.frame(tmpceRNA[ind1,],ratio[ind1])
#cesig_1 <- as.matrix(do.call(rbind,lapply(lapply(seq_len(dim(ceRNA_O_1)[1]), function(i) cesig[intersect(which(cesig[,1]==ceRNA_O_1[i,1]),which(cesig[,2]==ceRNA_O_1[i,2])),]),data.frame)))
#cenotsig_1 <- as.matrix(do.call(rbind,lapply(lapply(seq_len(dim(ceRNA_O_1)[1]), function(i) cenotsig[intersect(which(cenotsig[,1]==ceRNA_O_1[i,1]),which(cenotsig[,2]==ceRNA_O_1[i,2])),]),data.frame)))
 colnames(cesig)<-colnames(cenotsig)<-c("targetce","anotherce","miRNAs","miRNAs_num","ratio")
 rownames(cesig)<-rownames(cenotsig)<-NULL
}else if(method=="hypergeometric"){
allmiR_l<-length(unique(miRtar[,1]))
pvalue<-1-phyper(tmpl_shared-1,tmpl_miRs1,allmiR_l-tmpl_miRs1,tmpl_miRs2)
ind<-pvalue<cutoff
ind1<-!ind
cesig<-data.frame(tmpceRNA[ind,],pvalue[ind])
cenotsig<-data.frame(tmpceRNA[ind1,],pvalue[ind1])
#cesig_1 <- as.matrix(do.call(rbind,lapply(lapply(seq_len(dim(ceRNA_O_1)[1]), function(i) cesig[intersect(which(cesig[,1]==ceRNA_O_1[i,1]),which(cesig[,2]==ceRNA_O_1[i,2])),]),data.frame)))
#cenotsig_1 <- as.matrix(do.call(rbind,lapply(lapply(seq_len(dim(ceRNA_O_1)[1]), function(i) cenotsig[intersect(which(cenotsig[,1]==ceRNA_O_1[i,1]),which(cenotsig[,2]==ceRNA_O_1[i,2])),]),data.frame)))
 colnames(cesig)<-colnames(cenotsig)<-c("targetce","anotherce","miRNAs","miRNAs_num" ,"pvalue")
 rownames(cesig)<-rownames(cenotsig)<-NULL
}#for else
result<-list(cesig=cesig,cenotsig=cenotsig)
return(result)
}
