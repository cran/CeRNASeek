ceRNA.cmi <-
function(miRtar,targetce=NULL,geneexp,miRexp,numMIR=1,num_perm=100,cutoff=0.05){
         geneexp <- as.matrix(geneexp)
         miRexp <- as.matrix(miRexp)
		     miRtar <- miRtar[intersect(which(miRtar[, 1] %in% rownames(miRexp)),which(miRtar[, 2] %in% rownames(geneexp))), ]
         miRtar_d<-split(miRtar[,1],miRtar[,2])
         miRs_name<-names(miRtar_d)
         mi_cond<-function(data){
             data=as.matrix(data)
             geneexp_g1<-sort(data[,1],index.return=TRUE)
             mirexp_1<-sort(data[,2],index.return=TRUE)
             geneexp_g2<-sort(data[,3],index.return=TRUE)
             index<-cbind( as.data.frame(geneexp_g1[2]), as.data.frame(mirexp_1[2]), as.data.frame(geneexp_g2[2]))
             ydat<-index
             N1<-dim(data)[1]
             N2<-dim(data)[2]
             for(i in 1:N2){
                 ydat[index[,i],i]=matrix(1:N1,ncol=1);
                }
             dim1=2^(N2); dim2=2*(N2);
             rm(index);
             poc<-c();kon<-c();
             xcor=0; npar=1; poc[1]=1; kon[1]=N1; poradi=1:N1;
             NN<-matrix(0,nrow=1,ncol=dim1);marg=matrix(0,nrow=8*dim1,ncol=dim2);
             marg[1,]<-cbind(matrix(1,nrow=1,ncol=N2),matrix(N1,nrow=1,ncol=N2));
             Imm<-matrix(0:1,nrow=2,ncol=1);
             for(d in 2:N2){
                 Imm<-rbind(cbind(matrix(0,nrow=dim(Imm)[1],ncol=1),Imm),cbind(matrix(1,nrow=dim(Imm)[1],ncol=1),Imm));
                }
             chi2<-c(0,7.81,13.9,25.0,42.0);
             run=0;
             while(npar>0){
                 run=run+1;
                 apor<-c();
                 apoc<-poc[npar];akon<-kon[npar];
                 apor<-poradi[apoc:akon];Nex<-length(apor);
                 ave<-floor((marg[npar,1:N2]+marg[npar,(N2+1):dim2])/2);
                 J<-(ydat[apor,]<= matrix(1,ncol=1,nrow=Nex)%*%t(as.matrix(ave)))
                 I<-matrix(0,nrow=Nex,ncol=dim1)
                 amarg<-matrix(1,ncol=1,nrow=dim1)%*%marg[npar,];
                 for(d in 1:dim1){
                     I[,d]<-matrix(1,ncol=1,nrow=Nex)
                     for(k in 1:N2){
                         if(Imm[d,k]){
                             I[,d]<-as.numeric(I[,d]&!J[,k]);
                             amarg[d,k]<-ave[k]+1;
                            }
                         else{
                             I[,d]<-as.numeric(I[,d]&J[,k]);
                             amarg[d,k+N2]<-ave[k];
                            } 
                        }
                    }
                 NN<-apply(I,2,sum) 
                 tst<-dim1*sum((NN-Nex/dim1*(matrix(1,nrow=1,ncol=dim1)))^2)/Nex;
                 if((tst>chi2[N2])|run==1){
                     npar=npar-1;
                     for(ind in 1:dim1){
                         if(NN[ind]>dim1){
                             npar<-npar+1;
                             akon<-apoc+NN[ind]-1;
                             poc[npar]<-apoc;kon[npar]<-akon;
                             marg[npar,]<-amarg[ind,];
                             poradi[apoc:akon]<-apor[which(I[,ind]>0)];
                             apoc<-akon+1;
                            }
                         else{
                             if(NN[ind]>0){
                                 Nxx<-prod(amarg[ind,(N2+1):dim2]-amarg[ind,1:N2]+matrix(1,nrow=1,ncol=N2));
                                 Nz<-amarg[ind,6]-amarg[ind,3]+1;
                                 Jx <- as.numeric((ydat[,1]>=amarg[ind,1])&(ydat[,1]<=amarg[ind,4]));
                                 Jy <- as.numeric((ydat[,2]>=amarg[ind,2])&(ydat[,2]<=amarg[ind,5]));
                                 Jz <- as.numeric((ydat[,3]>=amarg[ind,3])&(ydat[,3]<=amarg[ind,6]));
                                 Nxz <- sum(as.numeric(Jx&Jz));
                                 Nyz <- sum(as.numeric(Jy&Jz));
                                 cond <- (NN[ind]*Nz)/(Nxz*Nyz);
                                 if(is.infinite(cond)){
                                     cond = 1;
                                    }
                                 else if(cond==0){
                                     cond = 1;
                                    }
                                 xcor = xcor + NN[ind]*log(cond);  
                                }
                            }
                        }
                    }
                 else{
                     Nxx <- prod(marg[npar,(N2+1):dim2]-marg[npar,1:N2]+matrix(1,nrow=1,ncol=N2));
                     Nz <- marg[npar,6]-marg[npar,3]+1;
                     Jx <- as.numeric((ydat[,1]>=marg[npar,1])&(ydat[,1]<=marg[npar,4]));
                     Jy <- as.numeric((ydat[,2]>=marg[npar,2])&(ydat[,2]<=marg[npar,5]));
                     Jz <- as.numeric((ydat[,3]>=marg[npar,3])&(ydat[,3]<=marg[npar,6]));
                     Nxz <- sum(as.numeric(Jx&Jz));
                     Nyz <- sum(as.numeric(Jy&Jz));
                     cond <- (Nex*Nz)/(Nxz*Nyz);
                     if(is.infinite(cond)){
                         cond = 1;
                        }
                     else if(cond==0){
                         cond = 1;
                        }
                     xcor = xcor + Nex*log(cond);
                     npar=npar-1;
                    }
                }
             xcor=xcor/N1;
             return(xcor);
            }
         if(is.null(targetce)){ 
             #ceRNA_o<-permutations(length(miRs_name),2,miRs_name,repeats.allowed=F)
			 #ceRNA_O_1<-combinations(length(miRs_name),2,miRs_name,repeats.allowed=F)
			 t(combn(miRs_name,2))-> ceRNA_O_1
             fun<-function(x){
                 tmp<-intersect(miRtar_d[[x[1]]],miRtar_d[[x[2]]])
                 tmp_p<-paste(tmp,collapse=",")
                 tmp_l<-length(tmp)
                 return(c(tmp_p,tmp_l))    
                }#for fun
             miRs<-t(apply(ceRNA_O_1,1,fun))
			 
             ceRNA_1<-data.frame(ceRNA_O_1,miRs[,1],as.numeric(miRs[,2]),stringsAsFactors=FALSE);
             ceRNA<-ceRNA_1[which(ceRNA_1[,4]!=0),]
			 #ceRNA_uni <- as.matrix(do.call(rbind,lapply(lapply(seq_len(dim(ceRNA_O_1)[1]), function(i) ceRNA[intersect(which(ceRNA[,1]==ceRNA_O_1[i,1]),which(ceRNA[,2]==ceRNA_O_1[i,2])),]),data.frame)))
			 ceRNA_F<-ceRNA[(as.numeric(ceRNA[,4]))>=numMIR,]
             fun1<-function(x){
                     sum=0
                     index_g1<-x[1]
                     index_g2<-x[2]
                     index_miR<-unlist(strsplit(x[3],split=","))
                     geneexp_g1<-geneexp[index_g1,];geneexp_g2<-geneexp[index_g2,];
					 Exp_data <- lapply(seq_len(length(index_miR)), function(i) cbind(geneexp_g1,miRexp[index_miR[i],],geneexp_g2))
                     CMI <-unlist( lapply(seq_len(length(Exp_data)), function(i) mi_cond(Exp_data[[i]])))
					 tmp <- as.matrix(do.call(rbind,lapply(lapply(seq_len(length(index_miR)), function(i) cbind(index_g1,index_miR[i],index_g2,CMI[[i]])),data.frame)))
                } 
             result<- as.matrix(do.call(rbind,lapply(apply(ceRNA_F,1,fun1), data.frame)));
             fun.random<-function(x){
                 index_g1<-x[1];index_g2<-x[3];index_miR<-x[2];
                 geneexp_g1<-geneexp[index_g1,];miRexp_m<-miRexp[index_miR,];
				 rand_Exp <- t(sapply(seq_len(num_perm), function(i) sample(geneexp[index_g2,])))
				 Exp_data <- lapply(seq_len(num_perm), function(i) cbind(geneexp_g1,miRexp_m,rand_Exp[i,]))
                 CMI <-unlist( lapply(seq_len(num_perm), function(i) mi_cond(Exp_data[[i]])))
				 tmp_rand<- as.matrix(do.call(rbind,lapply(lapply(seq_len(num_perm), function(i) t(CMI[[i]])),matrix)))
                 pvalue<-max(1,sum( as.numeric(x[4]<tmp_rand)))/num_perm ;
                 return(pvalue)
                }
			
				
             cmi_p<-t(apply(result,1,fun.random))
             ceRNA_cmi<-data.frame(result,t(cmi_p),stringsAsFactors=FALSE)
             rownames(ceRNA_cmi)<-NULL
             colnames(ceRNA_cmi)<-c("targetce","miRNAs","anotherce","cmi","pvalue")
             combP.fun <- function(data){
                 xsq<-(-2*sum(log(data)));comp<-c()
                 comp <- pchisq( xsq, 2 * length(data), lower.tail = FALSE);
                 return (c(xsq,comp))
                }
              p_list <- lapply(seq_len(dim(ceRNA_F)[1]), function(i) ceRNA_cmi[intersect(which(ceRNA_F[i,1]==ceRNA_cmi[,1]),which(ceRNA_F[i,2]==ceRNA_cmi[,3])),5])
              comP <- lapply(seq_len(dim(ceRNA_F)[1]), function(i) combP.fun(p_list[[i]]))
              ceRNA_comP<- cbind(ceRNA_F,do.call(rbind, as.matrix(comP)))
             rownames(ceRNA_comP)<-NULL
             colnames(ceRNA_comP)<-c("targetce","anotherce","commonMiR","number_miR","xsq","comP")
			 ceRNA_comP<-ceRNA_comP[ceRNA_comP[,6]<=cutoff,]
             result<-list(ceRNA_cmi=ceRNA_cmi,ceRNA_comP=ceRNA_comP);
             return(result)
            }
         else{
             tarmiRs<-miRtar_d[[targetce]]
             fun3<-function(x){
                 tmp<-intersect(miRtar_d[[x]],tarmiRs)
                 tmp_p<-paste(tmp,collapse=",")
                 tmp_l<-length(tmp)
                 return(c(tmp_p,tmp_l))
                }#for fun1
             miRs<-t(sapply(miRs_name,fun3))
             ceRNA<-data.frame(targetce,rownames(miRs),miRs[,1],as.numeric(miRs[,2]),stringsAsFactors=FALSE)
             ceRNA<-ceRNA[intersect(which(ceRNA[,4]!=0),which(ceRNA[,1]!=ceRNA[,2])),]
			 ceRNA_uni<-ceRNA[(as.numeric(ceRNA[,4]))>=numMIR,]
             rm(miRs);
             fun4<-function(x){
                     sum=0
                     index_g1<-x[1]
                     index_g2<-x[2]
                     index_miR<-unlist(strsplit(x[3],split=","))
                     geneexp_g1<-geneexp[index_g1,];geneexp_g2<-geneexp[index_g2,];
					 Exp_data <- lapply(seq_len(length(index_miR)), function(i) cbind(geneexp_g1,miRexp[index_miR[i],],geneexp_g2))
                     CMI <-unlist( lapply(seq_len(length(Exp_data)), function(i) mi_cond(Exp_data[[i]])))
					 tmp <- as.matrix(do.call(rbind,lapply(lapply(seq_len(length(index_miR)), function(i) cbind(index_g1,index_miR[i],index_g2,CMI[[i]])),data.frame)))
                } 
             result<- as.matrix(do.call(rbind,lapply(apply(ceRNA_uni,1,fun4), data.frame)));
            fun.random<-function(x){
                 index_g1<-x[1];index_g2<-x[3];index_miR<-x[2];
                 geneexp_g1<-geneexp[index_g1,];miRexp_m<-miRexp[index_miR,];
				 rand_Exp <- t(sapply(seq_len(num_perm), function(i) sample(geneexp[index_g2,])))
				 Exp_data <- lapply(seq_len(num_perm), function(i) cbind(geneexp_g1,miRexp_m,rand_Exp[i,]))
                 CMI <-unlist( lapply(seq_len(num_perm), function(i) mi_cond(Exp_data[[i]])))
				 tmp_rand<- as.matrix(do.call(rbind,lapply(lapply(seq_len(num_perm), function(i) t(CMI[[i]])),matrix)))
                 pvalue<-max(1,sum( as.numeric(x[4]<tmp_rand)))/num_perm ;
                 return(pvalue)
                }
             cmi_p<-t(apply(result,1,fun.random))
             ceRNA_cmi<-data.frame(result,t(cmi_p),stringsAsFactors=FALSE)
             rownames(ceRNA_cmi)<-NULL
             colnames(ceRNA_cmi)<-c("targetce","miRNA","anotherce","cmi","pvalue") 
			 
			 
			 combP_fun <- function(data){
             xsq<-(-2*sum(log(data)));comp<-c()
             comp <- pchisq( xsq, 2 * length(data), lower.tail = FALSE);
             return (c(xsq,comp))
              }
              p_list <- lapply(seq_len(dim(ceRNA_uni)[1]), function(i) ceRNA_cmi[intersect(which(ceRNA_uni[i,1]==ceRNA_cmi[,1]),which(ceRNA_uni[i,2]==ceRNA_cmi[,3])),5])
              comP <- lapply(seq_len(dim(ceRNA_uni)[1]), function(i) combP_fun(p_list[[i]]))
              ceRNA_comP<- cbind(ceRNA_uni,do.call(rbind, as.matrix(comP)))
			 rownames(ceRNA_comP)<-NULL
             colnames(ceRNA_comP)<-c("targetce","anotherce","commonMiR","number_miR","xsq","comP") 
            }
		  ceRNA_comP<-ceRNA_comP[ceRNA_comP[,6]<=cutoff,]#for else
          result<-list(ceRNA_cmi=ceRNA_cmi,ceRNA_comP=ceRNA_comP);
          return(result)
    }
