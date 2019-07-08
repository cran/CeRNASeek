ceRNA.enrich <-
function(data,GOterms,background,threshold=2,correction="BH"){
  #filter the GOterms genes to make sure that all the Goterms genes are in the background
    goname<-names(GOterms)
    index<-numeric()
    for(i in 1:length(goname)){
      comm<-intersect(GOterms[[goname[i]]],background)
      if(length(comm)!=0){
        assign(goname[i],comm)
        index<-c(index,i)
      }
    }
    goname<-goname[index]
  
  #an object named "goterm_freq" was constructed to contain the gene number in goterms
  goterm_num<-as.data.frame(matrix(nrow=length(goname),ncol=2))
  colnames(goterm_num)<-c("goterm","go_num")
  goterm_num[,1]<-goname
  for(j in 1:length(goname)){goterm_num[j,2]<-length(get(goname[j]))}
  
  #an object named "lnc_freq" was constructed to contain the target gene number of an lncRNA(or more lncRNAs)
  #predict the goterms enriched by the target genes of an lncRNA(or more lncRNAs)
  tar<-sort(as.character(unique(data[,1])))
  tar_num<-as.data.frame(matrix(ncol=2,nrow=length(tar)))
  tar_num[,1]<-tar
  colnames(tar_num)<-c("tar","tarnum")
  inter<-list()
  for(j in 1:length(tar)){
    tmptar<-unique(data[data[,1]==tar[j],2])
    tmptar_num<-length(tmptar)
    tmpinter<-numeric()
    for(k in goname){
      tmpinter<-append(tmpinter,length(intersect(tmptar,get(k))))
    }#for k(GOterm name)
    inter[[j]]<-tmpinter
    tar_num[j,2]<-length(tmptar)
  }#for j(the index of lncRNA)
  interd<-unlist(inter)
  tard<-sort(rep(tar,length(goname)))
  gotermd<-rep(goname,length(tar))
  result<-as.data.frame(cbind(tard,gotermd))
  result<-cbind(result,interd)
  colnames(result)<-c("tar","goterm","internum")
  result<-result[result[,"internum"]>=threshold,]
  result<-merge(result,goterm_num)
  result<-merge(result,tar_num)
  
  #p-value adjusted by the method assigned
  backnum<-length(background)
  pvalue<-apply(result[,c("tarnum","go_num","internum")],1,function(x){return(1-phyper(x[3],x[1],backnum-x[1],x[2]))})
  result<-result[,c("tar","goterm","tarnum","go_num","internum")]
  result<-cbind(result,pvalue)
  fdr<-p.adjust(result[,"pvalue"],method=correction)
  result<-cbind(result,fdr)
  colnames(result)<-c("target","GOterm","target_num","GOtermnum","term_tar","P_value","fdr")
  return(result)
}
