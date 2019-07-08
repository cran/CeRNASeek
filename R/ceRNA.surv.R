ceRNA.surv <-
function(ceRNA,exp.sur,train=NULL,test=NULL,index){
 if(is.null(train)&is.null(train)){
     wr<-as.data.frame(matrix(ncol=10,nrow=length(index)),stringsAsFactors=F)
     colnames(wr)<-c("targetce","anotherce","coef_targetce","p_targetce","coef_anotherce","p_anotherce","N_low","N_high","HR","p")
    for(i in index){
  #calculation and plot
    exp_sur<-matrix(nrow=dim(exp.sur)[1],ncol=2) 
    exp_sur[,1]<-exp.sur[,ceRNA[i,1]]
    exp_sur[,2]<-exp.sur[,ceRNA[i,2]]
    colnames(exp_sur)<-ceRNA[i,1:2]
    exp_sur<-t(exp_sur)
	
    q<-summary(coxph(Surv(OS,living)~exp.sur[,ceRNA[i,1]],exp.sur))
    wr[,1]=ceRNA[i,1]
    wr[,3]=q$coefficients[,1] 
    wr[,5]=q$coefficients[,5]
    w1<-q$coefficients[,2]	
    q<-summary(coxph(Surv(OS,living)~exp.sur[,ceRNA[i,2]],exp.sur))
    wr[,2]=ceRNA[i,2]
    wr[,4]=q$coefficients[,1] 
    wr[,6]=q$coefficients[,5] 
	w2<-q$coefficients[,2]	
  
    
    exp_score<-as.data.frame(matrix(ncol=4,nrow=dim(exp.sur)[1]))
    exp_score[,1]<-rownames(exp_sur)
    exp_score[,2]<-exp.sur[,"living"]
    exp_score[,3]<-exp.sur[,"OS"]
    exp_score[,4]<-(w1/(w1+w2))*exp_sur[1,]+(w2/(w1+w2))*exp_sur[2,]
    cuff<-sort(exp_score[,4],decreasing=F)
    index_exp<-order(exp_score[,4])[-(ceiling(length(cuff)*0.25):(length(cuff)-(ceiling(length(cuff)*0.25))+1))]
    #cutoff<-train_score[ceiling(length(cuff_train<-sort(train_score[,4],decreasing=F))*0.25),]
    exp_top<-as.data.frame(matrix(ncol=5,nrow=length(index_exp)))
    exp_top[,1]<-exp_score[index_exp,1]
    exp_top[,2]<-exp_score[index_exp,2]
    exp_top[,3]<-exp_score[index_exp,3]
    exp_top[,4]<-exp_score[index_exp,4]
    
    index1<-exp_top[,4]<=sort(exp_score[,4],decreasing=F)[ceiling(length(cuff)*0.25)-1]
    index2<-!(index1)
    exp_top[index1,5]<-1
    exp_top[index2,5]<-2
    
    index1_sum<-sum(index1)
    index2_sum<-sum(index2)
    wr[,7]<-index1_sum
    wr[,8]<-index2_sum
    dif<-survdiff(Surv(exp_top[,3],exp_top[,2])~exp_top[,5],exp_top)
    p <- 1-pchisq(dif$chisq,length(dif$n) - 1)
	HR <- (dif$obs[1]/dif$exp[1])/(dif$obs[2]/dif$exp[2])
	wr[,9]<-HR
    wr[,10]<-p
	surfit<-survfit(Surv(exp_top[,3],exp_top[,2])~exp_top[,5],exp_top)
	title<-paste("exp:",paste(ceRNA[i,1:2],collapse = "-"))
    plot(surfit,col=c('#5D83AA','#C4312F'),lwd=2,xlab="Survival time(days)",ylab="Survival probability",main=title)
    p<-paste("p-value=",round(p,3))
    text(max(exp_score[,3],na.rm = T)*0.5,0.9,p,pos=4)
    legend("topright",c("LowRisk","HighRisk"),lty=1:1,col=c('#5D83AA','#C4312F'),bty = "n")
    
   }
   return(wr)
   
 }else{
  #a dataframe named "wr" was constructed to contain result
  wr<-as.data.frame(matrix(ncol=14,nrow=length(index)),stringsAsFactors=F)
  colnames(wr)<-c("targetce","anotherce","coef_targetce","p_targetce","coef_anotherce","p_anotherce","N_train_low","N_train_high","HR_train","p_train","N_test_low","N_test_high","HR_test","p_test")
   for(i in index){
  #calculation and plot
    exp_train<-matrix(nrow=dim(train)[1],ncol=2) 
    train1<-as.vector(t(train))	
    exp_train[,1]<-exp.sur[train1,ceRNA[i,1]]
    exp_train[,2]<-exp.sur[train1,ceRNA[i,2]]
    colnames(exp_train)<-ceRNA[i,1:2]
    rownames(exp_train)<-train1
    exp_train<-t(exp_train)
    exp_test<-matrix(nrow=dim(test)[1],ncol=2)
	test1<-as.vector(t(test))
    exp_test[,1]<-exp.sur[test1,ceRNA[i,1]]
    exp_test[,2]<-exp.sur[test1,ceRNA[i,2]]
    colnames(exp_test)<-ceRNA[i,1:2]
    rownames(exp_test)<-test1
    exp_test<-t(exp_test)
    #
    q<-summary(coxph(Surv(OS,living)~exp.sur[train1,ceRNA[i,1]],exp.sur[train1,]))
    wr[,1]=ceRNA[i,1]
    wr[,3]=q$coefficients[,1] 
    wr[,5]=q$coefficients[,5]
    w1<-q$coefficients[,2]	
    q<-summary(coxph(Surv(OS,living)~exp.sur[train1,ceRNA[i,2]],exp.sur[train1,]))
    wr[,2]=ceRNA[i,2]
    wr[,4]=q$coefficients[,1] 
    wr[,6]=q$coefficients[,5] 
	w2<-q$coefficients[,2]	
  
    
    train_score<-as.data.frame(matrix(ncol=4,nrow=dim(train)[1]))
    train_score[,1]<-train
    train_score[,2]<-exp.sur[train1,"living"]
    train_score[,3]<-exp.sur[train1,"OS"]
    train_score[,4]<-(w1/(w1+w2))*exp_train[1,]+(w2/(w1+w2))*exp_train[2,]
    cuff_train<-sort(train_score[,4],decreasing=F)
    index_train<-order(train_score[,4])[-(ceiling(length(cuff_train)*0.25):(length(cuff_train)-(ceiling(length(cuff_train)*0.25))+1))]
    #cutoff<-train_score[ceiling(length(cuff_train<-sort(train_score[,4],decreasing=F))*0.25),]
    train_top<-as.data.frame(matrix(ncol=5,nrow=length(index_train)))
    train_top[,1]<-train_score[index_train,1]
    train_top[,2]<-train_score[index_train,2]
    train_top[,3]<-train_score[index_train,3]
    train_top[,4]<-train_score[index_train,4]
    
    index1<-train_top[,4]<=sort(train_score[,4],decreasing=F)[ceiling(length(cuff_train)*0.25)-1]
    index2<-!(index1)
    train_top[index1,5]<-1
    train_top[index2,5]<-2
    
    index1_sum<-sum(index1)
    index2_sum<-sum(index2)
    wr[,7]<-index1_sum
    wr[,8]<-index2_sum
   
     par(mfrow=c(length(index),2))
    dif<-survdiff(Surv(train_top[,3],train_top[,2])~train_top[,5],train_top)
    p <- 1-pchisq(dif$chisq,length(dif$n) - 1)
	HR <- (dif$obs[1]/dif$exp[1])/(dif$obs[2]/dif$exp[2])
	wr[,9]<-HR
    wr[,10]<-p
	surfit<-survfit(Surv(train_top[,3],train_top[,2])~train_top[,5],train_top)
	title<-paste("train:",paste(ceRNA[i,1:2],collapse = "-"))
    plot(surfit,col=c('#5D83AA','#C4312F'),lwd=2,xlab="Survival time(days)",ylab="Survival probability",main=title)
    p<-paste("p-value=",round(p,3))
    text(max(train_score[,3],na.rm = T)*0.5,0.9,p,pos=4)
    legend("topright",c("LowRisk","HighRisk"),lty=1:1,col=c('#5D83AA','#C4312F'),bty = "n")
    
    test_score<-as.data.frame(matrix(ncol=4,nrow=dim(test)[1]))
    test_score[,1]<-test
    test_score[,2]<-exp.sur[test1,"living"]
    test_score[,3]<-exp.sur[test1,"OS"]
    test_score[,4]<-(w1/(w1+w2))*exp_test[1,]+(w2/(w1+w2))*exp_test[2,]
    cuff_train<-sort(test_score[,4],decreasing=F)
	index_test<-order(test_score[,4])[-(ceiling(length(cuff_train)*0.25):(length(cuff_train)-(ceiling(length(cuff_train)*0.25))+1))]
	#cutoff<-test_score[ceiling(length(cuff_train<-sort(test_score[,4],decreasing=F))*0.25),]
	
	test_top<-as.data.frame(matrix(ncol=5,nrow=length(index_test)))
	test_top[,1]<-test_score[index_test,1]
    test_top[,2]<-test_score[index_test,2]
    test_top[,3]<-test_score[index_test,3]
    test_top[,4]<-test_score[index_test,4]
    
    index1<-test_top[,4]<=sort(train_score[,4],decreasing=F)[ceiling(length(cuff_train)*0.25)-1]
    index2<-!(index1)
    test_top[index1,5]<-1
    test_top[index2,5]<-2
    
    index1_sum<-sum(index1)
    index2_sum<-sum(index2)
    wr[,11]<-index1_sum
    wr[,12]<-index2_sum
    
    dif<-survdiff(Surv(test_top[,3],test_top[,2])~test_top[,5],test_top)
	HR <- (dif$obs[1]/dif$exp[1])/(dif$obs[2]/dif$exp[2])
    p <- 1-pchisq(dif$chisq,length(dif$n) - 1)
	wr[,13]<-HR
    wr[,14]<-p
    survfit<-survfit(Surv(test_top[,3],test_top[,2])~test_top[,5],test_top)
    title<-paste("test:",paste(ceRNA[i,1:2],collapse = "-"))
    plot(surfit,col=c('#5D83AA','#C4312F'),lwd=2,xlab="Survival time(days)",ylab="Survival probability",main=title)
    p<-paste("p-value=",round(p,3))
    text(max(train_score[,3],na.rm = T)*0.5,0.9,p,pos=4)
    legend("topright",c("LowRisk","HighRisk"),lty=1:1,col=c('#5D83AA','#C4312F'),bty = "n")

    
   }
   return(wr)
  }
}
