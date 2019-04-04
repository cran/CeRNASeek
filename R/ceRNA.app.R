ceRNA.app <-
function(miRtar,targetce=NULL){
miRtar_d<-split(miRtar[,1],miRtar[,2])
miR_l<-sapply(miRtar_d,length)

miRs_name<-names(miRtar_d)
if(is.null(targetce)){
ceRNA_o<-permutations(length(miRs_name),2,miRs_name,repeats.allowed=T)######permutations
fun<-function(x){
tmp<-intersect(miRtar_d[[x[1]]],miRtar_d[[x[2]]])
tmp_p<-paste(tmp,collapse=",")
tmp_l<-length(tmp)
return(c(tmp_p,tmp_l))
}#for fun
miRs<-t(apply(ceRNA_o,1,fun))
ceRNA<-data.frame(ceRNA_o,miRs[,1],as.numeric(miRs[,2]),stringsAsFactors=FALSE)
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
