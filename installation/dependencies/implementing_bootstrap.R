###### LOADING LIBRARIES #####
# Load libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(dendextend)
library(cluster)


###### ADAPTED FUNCTIONS #####
# adapted from the package "bactaxR" (paper: doi 10.1128/mBio.00034-20)
#  
# read.ANI.custom custom function
read.ANI.custom <- function (file) {
  f <- read.delim(file = file, header = F, sep = "\t", stringsAsFactors = F)
  f <- f[, 1:3]
  colnames(f) <- c("query", "reference", "ANI")
  f$query <- unlist(lapply(strsplit(x = f$query, split = "/"), 
                           function(x) x[length(x)]))
  f$reference <- unlist(lapply(strsplit(x = f$reference, split = "/"), 
                               function(x) x[length(x)]))
  if (all(f$query %in% f$reference)) {
    return(f)
  }
  else {
    stop("Error: file does not contain pairwise comparisons for all genomes.")
  }
}

# ANI.dendrogram.custom custom function
ANI.dendrogram.custom <- function (ani_object, ANI_threshold = 95, label_size = 1e-10, xline = NULL, xlinecol = "#20A387FF", 
                                   xlinetype = "dashed") {
  
  ani_object <- ani
  ANI_threshold = 95
  label_size = 1e-10
  xline = NULL
  xlinecol = "#20A387FF"
  xlinetype = "dashed"
  
  fastani <- ani_object
  colnames(fastani) <- c("query", "reference", "ANI")
  s <- dcast(fastani, formula <- query ~ reference, value.var = "ANI")
  rownames(s) <- s$query
  s <- as.matrix(s[, !(colnames(s) == "query")])
  j <- matrix(data = 100, nrow = nrow(s), ncol = ncol(s))
  d <- j - s
  d <- 0.5 * (d + t(d))
  #d <- as.dist(d)
  
  dism <- function(x){
    res <- as.dist(x)
    print(dim(x))
    attr(res, "method") <- "dism"
    return(res)
  }
  
  #d.sym <- 0.5 * (d + t(d))
  #d.dist <- as.dist(d.sym)
  #h <- hclust(d = d.dist, method = "average")
  #h <- pvclust(d = d.sym, method.dist="cor", method.hclust="average", nboot=10)
  h <- pvclust(d = d, method.dist=dism, method.hclust="average", nboot=10)
  
  dend <- as.dendrogram(h)
  par(mar = c(10,10,10,10))
  #dend <- dend %>% set("labels_cex", label_size) 
  dend <- dend %>% set("labels_cex", label_size)
  
  plot(dend, horiz = T, xlab = "Dissimilarity (100-ANI)", main = taxonomy, xlim=c(20,0)) %>% text

  
  if (!(is.null(xline))) {
    abline(v = xline, lty = xlinetype, col = xlinecol)
  }
  if (!(is.null(ANI_threshold))) {
    clusters <- cutree(tree = dend, h = 100 - ANI_threshold)
    nclus <- length(table(clusters))
    medoid.genome <- c()
    medoid.cluster <- c()
    allgenomes.cluster <- c()
    allgenomes.genome <- c()
    allgenomes.medoid <- c()
    for (i in 1:nclus) {
      genomes <- clusters[which(clusters == i)]
      if (length(genomes) > 1) {
        genomes.d <- d.sym[which(rownames(d.sym) %in% 
                                   names(genomes)), ]
        genomes.d <- genomes.d[, which(colnames(genomes.d) %in% 
                                         names(genomes))]
        gnames <- rownames(genomes.d)
        genomes.d <- as.dist(genomes.d)
        medoid <- pam(x = genomes.d, k = 1)
        medoid.genome <- c(medoid.genome, medoid$medoids)
        allgenomes.genome <- c(allgenomes.genome, gnames)
        allgenomes.cluster <- c(allgenomes.cluster, rep(i, 
                                                        length(gnames)))
        allgenomes.medoid <- c(allgenomes.medoid, rep(medoid$medoids, 
                                                      length(gnames)))
      }
      else {
        medoid <- names(genomes)
        medoid.genome <- c(medoid.genome, medoid)
        allgenomes.cluster <- c(allgenomes.cluster, i)
        allgenomes.genome <- c(allgenomes.genome, medoid)
        allgenomes.medoid <- c(allgenomes.medoid, medoid)
      }
      medoid.cluster <- c(medoid.cluster, i)
    }
    medoid_genomes <- data.frame(medoid.cluster, medoid.genome)
    colnames(medoid_genomes) <- c("Cluster", "Genome")
    cluster_assignments <- data.frame(allgenomes.cluster, 
                                      allgenomes.genome, allgenomes.medoid)
    colnames(cluster_assignments) <- c("Cluster", "Genome", 
                                       "Cluster_Medoid")
    
    clusters_order <- cutree(tree = dend, h = 100 - ANI_threshold, order_clusters_as_data = FALSE)
    clusters_order <- as.data.frame(clusters_order)
    
    dend %>% color_branches(dend, k=nclus, groupLabels = T) %>% color_labels(dend, k=nclus) %>%
      plot(horiz = T, xlab = "Dissimilarity (100-ANI)", main = taxonomy, xlim=c(20,0))
    
    if (!(is.null(xline))) {
      abline(v = xline, lty = xlinetype, col = xlinecol)
    }
    return(list(medoid_genomes = medoid_genomes, cluster_assignments = cluster_assignments))
  }
  else {
    return(NULL)
  }
}

plot_dendrogram <- function(){
  # Generate plot and return clusters
  cluster_dendrogram <- ANI.dendrogram.custom(ani_object = ani,
                                              ANI_threshold = 95,
                                              #xline = c(4,5,6,7.5),
                                              #xlinecol = c("#ffc425", "#f37735", "deeppink4", "black"),
                                              xline = 5,
                                              xlinecol = "deeppink4",
                                              label_size = 0.5)
}

#setwd("Desktop/tests_dendro/cluster-95/")
#ani <- read.ANI.custom("fastani-out-1500-0.txt")
#taxonomy <- "test"

plot_dendrogram()
