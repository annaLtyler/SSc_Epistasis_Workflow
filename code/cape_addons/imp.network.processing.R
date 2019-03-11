imp.network.processing <- function(){

library(igraph)
full.net.file <- "source_target_net_chr4_chr18_expanded.csv"
source.gene.file <- "chr4_genes.txt"
target.gene.file <- "chr18_genes.txt"

edge.list <- as.matrix(read.table(full.net.file, header = TRUE, stringsAsFactors = FALSE, sep = ","))

net <- graph.edgelist(matrix(edge.list[,1:2], ncol = 2), directed = FALSE)
E(net)$weight <- as.numeric(edge.list[,3])

#make sure you capitalize everything in these files
source.nodes <- as.matrix(read.table(source.gene.file, stringsAsFactors = FALSE, sep = ","))
if(dim(source.nodes)[2] > 1){
	source.nodes <- matrix(unique(c(source.nodes[,1], source.nodes[,2])), ncol = 1)
	}
target.nodes <- as.matrix(read.table(target.gene.file, stringsAsFactors = FALSE, sep = ","))
if(dim(target.nodes)[2] > 1){
	target.nodes <- matrix(unique(c(target.nodes[,1], target.nodes[,2])), ncol = 1)
	}

source.locale <- which(V(net)$name %in% source.nodes[,1])
target.locale <- which(V(net)$name %in% target.nodes[,1])

node.colors <- rep(NA, vcount(net))

node.colors[source.locale] <- "red"
node.colors[target.locale] <- "blue"

V(net)$color <- node.colors

threshold <- seq(0.5, 0.99, 0.01)

thresh.net <- net
pdf("net.pdf")
for(th in 1:length(threshold)){
	print(th)
	low.edges <- which(E(thresh.net)$weight < threshold[th])
	thresh.net <- delete.edges(thresh.net, low.edges) 
	plot(thresh.net, vertex.size = 2, vertex.label = NA, main = paste("Threshold:", threshold[th], sep = " "))
	legend("topleft", legend = c("source", "target"), col = c("red", "blue"), pch = 16)
	}
dev.off()


choice.threshold = 0.93
num.clust = 1
low.edges <- which(E(net)$weight < choice.threshold)
thresh.net <- delete.edges(net, low.edges)

names.of.interest <- c(source.nodes[,1], target.nodes[,1])
source.target.names <- rep(NA, vcount(thresh.net))
names.order <- match(V(thresh.net)$name, names.of.interest)
names.locale <- which(V(thresh.net)$name %in% names.of.interest)
source.target.names[names.locale] <- names.of.interest[names.order[!is.na(names.order)]]

# net.layout <- NULL
# net.layout <- as.matrix(read.table("net.layout.txt"))
net.layout <- get.genes.in.cluster(thresh.net, num.clusters = num.clust, plot.clusters = TRUE, source.target.names, net.layout = net.layout)
net.layout <- get.genes.in.cluster(thresh.net, num.clusters = num.clust, plot.clusters = TRUE, v.labels = V(net)$name, net.layout = net.layout)
# write.table(net.layout, file = "net.layout.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

thresh.edges <- get.edgelist(thresh.net)
source.locale <- c(which(thresh.edges[,1] %in% source.nodes[,1]), which(thresh.edges[,2] %in% source.nodes[,1]))
target.locale <- c(which(thresh.edges[,1] %in% target.nodes[,1]), which(thresh.edges[,2] %in% target.nodes[,1]))
direct.connect <- thresh.edges[intersect(target.locale, source.locale),]
direct.connect


net.clust <- clusters(thresh.net)

max.clust.nodes <- V(thresh.net)$name[which(net.clust$membership == 1)]
write.table(max.clust.nodes, paste("Major_Component_Nodes_", choice.threshold, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
cat(length(max.clust.nodes), "genes in large component")


}