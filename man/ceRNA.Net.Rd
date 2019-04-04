\name{ceRNA.Net}
\alias{ceRNA.Net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visualize and analyze the identified ceRNA-ceRNA network using ceRNA.Net function 
}
\description{
A downstream analysis function of ceRNAseekvisualize and analyze the identified ceRNA-ceRNA network, the network can
 be defined as weighted or un-weighted network.and the basic topological features (such as degree, 
 closeness, betweenness and centrality) of the ceRNAs can be output.
}
\usage{
ceRNA.Net(data, net_direct = TRUE, vertex_size = 20, v.label = TRUE,
         node_shape = "circle",n_color = "orange",E_weight = TRUE, ity = 1,
         label_cex = 2, label_color = "black",edge_color = "gray", 
         n_frame.color = "gray")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix of ceRNA interaction pairs identified by statistical identification methods selected by the user .
}
  \item{net_direct}{
A logical variable specifies a directed or undirected network.default (net_direct = TRUE).
}
  \item{vertex_size}{
a numeric vector to adjust the node size,default (vertex_size = 20).
}
  \item{v.label}{
Whether to display the label of the node.default (v.label = TRUE).
}
  \item{node_shape}{
A character vector is used to adjust the shape of the node,The selectable shape parameters are
"circle","square","csquare","rectangle","crectangle","vrectangle","pie",
"sphere","none".default (node_shape = "circle").
}
  \item{n_color}{
The character vector used to define the fill color of the node.
}
  \item{E_weight}{
A logical vector represents whether the network is a weighted network,default (E_weight = TRUE).
}
  \item{ity}{
A numeric vector defines whether the edge is a solid line or a dotted line,the possible values of the vector are c(1, 2),
default (ity = 1).
}
  \item{label_cex}{
Specify the size of the node label font.
}
  \item{label_color}{
Specify the label color of the node.
}
  \item{edge_color}{
Specify the color of the network side.
}
  \item{n_frame.color}{
The character vector used to define the border color of the node.
}
}
\details{
This function calls the igraph package. For specific parameter settings, please refer to the igrap help documentation.
Note:All the arguments without default value must be assigned.
}
\value{
The output includes two parts:
the network diagram of ceRNA interaction and the topology attribute information of the network.
\itemize{
Network topology attributes include 5 types of information:
\itemize{
\item \code{degree} degree refers to the number of edges in the network directly connected to the node
\item \code{closeness} An indicator that describes the average distance of a node to all other nodes in the network.
\item \code{betweenness} The proportion of this node that appears in the shortest path between other nodes.
\item \code{cluster coefficient} Representing the dense connection nature between some nodes in the network
\item \code{Eigenvector centrality} Representing the characteristic vector centrality of the network.
}
}
}
\examples{
##Display ceRNA interactions in a network format and output network topology attributes.
##The input file can be a list of [ceRNA] of the ceRNA.Lin result file, a list of [cesig] 
##for the result file identified by ceRNA.basic, or a list of [ceRNA_comP] in the result
##file identified by the ceRNA.cmi function.
##Here we apply the ceRNA list in the example file for CMI identification to 
##display the network and analyze the network topology properties.
ceRNA.Net(as.matrix(dataset[["Pre.ceRNA"]]),net_direct=TRUE,vertex_size=20,v.label = TRUE,
         node_shape="circle",n_color = "orange",E_weight=TRUE,ity=1,label_cex=2,
         label_color="black",edge_color="gray",n_frame.color = "gray")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ceRNA.Net}
