ceRNA.Net <-
function(data,net_direct=TRUE,vertex_size=20,v.label=TRUE,node_shape="circle",n_color="orange",E_weight=TRUE,ity=1,label_cex=2,label_color="black",edge_color="gray",n_frame.color="gray"){
g <- graph.edgelist(data[,1:2],directed = net_direct)
 V(g)$label.cex <-label_cex; V(g)$label.color <- label_color;E(g)$color<-edge_color
 E(g)$lty<-ity;V(g)$color<-n_color; V(g)$frame.color <- n_frame.color
 tmp2<-closeness(g,mode="all")
 tmp1<-degree(g);tmp3<-betweenness(g);tmp4<-transitivity(g, type="local");tmp5<-evcent(g)$vector;
 ceRNA_Net<-data.frame(t(t(tmp1)),t(t(tmp2)),t(t(tmp3)),t(t(tmp4)),t(t(tmp5)))
 if(E_weight){
  fun.w<-function(x){
  miRs<-unlist(strsplit(x[3],split=","))
  len<-length(miRs)
  return(len)
  }
 weight<-t(apply(data,1,fun.w))
 #pdf(file=paste("ceRNA_weighted_interaction_network",".pdf",sep=""))
 E(g)$weight<-weight
 if(v.label){
 plot(g,vertex.size= vertex_size,vertex.label=V(g)$name,shape=node_shape,edge.width=E(g)$weight)
 }else{
 plot(g,vertex.size= vertex_size,vertex.label=NA,shape=node_shape,edge.width=E(g)$weight)
 }
# dev.off()
 }else{
# pdf(file=paste("ceRNA_unweighted_interaction_network",".pdf",sep=""))
 if(v.label){
 plot(g,vertex.size= vertex_size,vertex.label=V(g)$name,shape=node_shape)
 }else{
 plot(g,vertex.size= vertex_size,vertex.label=NA,shape=node_shape)
 }
# dev.off()
  }
 #cat("The plot(s) is on the working directory!\n\n")
 colnames(ceRNA_Net)<-c("degree","closeness","betweenness","cluster coefficient","Eigenvector centrality") 
 return(ceRNA_Net);
 }
