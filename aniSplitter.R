#!/usr/bin/env Rscript

library(optparse,lib.loc="/data/msb/tools/anisplitter/packages_new/")


###### ARGUMENT PARSING ######
suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
library("optparse")
option_list = list(
  make_option(c("-d", "--directory"),
              type="character",
              default=NULL, 
              help="Directory path",
              metavar="character"),
  make_option(c("-f", "--fastani"),
              type="character",
              default="fastani-out.txt",
              help="Fastani output file name [default=%default]",
              metavar="character"),
  make_option(c("-t", "--taxonomy"),
              type="character",
              default="cluster_taxonomy_gtdb",
              help="Taxonomy file name [default=%default]",
              metavar="character"),
  make_option(c("-i", "--iterations"),
              type="numeric",
              default=1000,
              help="Number of iterations for bootstrapping [default=%default]",
              metavar="character"),
  make_option(c("-s", "--seed"),
              type="numeric",
              default=NULL,
              help="Set seed for reproducibility [default=%default]",
              metavar="character"),
  make_option(c("-a", "--ani-threshold"),
              type="numeric",
              default=95, 
              help="ANI threshold for splitting into groups [default=%default]",
              metavar="numeric")
);

opt_parser = OptionParser(
  usage = "usage: %prog [options]",
  option_list=option_list,
  prog = "aniSplitter.R"
);

opt = parse_args(opt_parser);

if (is.null(opt$directory)){
  print_help(opt_parser)
  stop("At least one argument must be supplied", call.=FALSE)
}

# Define variables from parsed arguments
 directory <- opt$d # directory path
 fastani_output_file <- opt$f # fastani output file
 taxonomy_file <- opt$t # taxonomy file
 ani_threshold <- opt$a # ani threshold
 bootstrap_iterations <- opt$i # number of iterations for bootstrap
 seed <- opt$s # set.seed with this number

###### LOADING LIBRARIES #####
# Load libraries
message("Loading libraries ...")
#suppressPackageStartupMessages(library(reshape2),lib.loc="/data/msb/tools/anisplitter/packages/")
library(reshape2,lib.loc="/data/msb/tools/anisplitter/packages_new/")
#suppressPackageStartupMessages(library(ggplot2))
library(ggplot2,lib.loc="/data/msb/tools/anisplitter/packages_new/")
suppressPackageStartupMessages(library(dplyr))
#library(dplyr,lib.loc="/data/msb/tools/anisplitter/packages/")
#suppressPackageStartupMessages(library(dendextend))
library(dendextend,lib.loc="/data/msb/tools/anisplitter/packages_new/")
#suppressPackageStartupMessages(library(cluster))
library(cluster,lib.loc="/data/msb/tools/anisplitter/packages_new/")
#suppressPackageStartupMessages(library(fpc)) # clusterboot function
library(fpc,lib.loc="/data/msb/tools/anisplitter/packages_new/") # clusterboot function




 # functions adapted from the package "bactaxR" (paper: doi 10.1128/mBio.00034-20)
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
message("step 2")
 ANI.dendrogram.custom <- function (ani_object, ANI_threshold = 95,
                                    label_size = 1e-10) {
   
   # Load ANI object
   fastani <- ani_object
   colnames(fastani) <- c("query", "reference", "ANI")
   
   ####### Calculate ANI distances #######
   s <- dcast(fastani, formula <- query ~ reference, value.var = "ANI")
   rownames(s) <- s$query
   s <- as.matrix(s[, !(colnames(s) == "query")])
   j <- matrix(data = 100, nrow = nrow(s), ncol = ncol(s))
   d <- j - s
   # Make matrix symmetric
   d.sym <- 0.5 * (d + t(d))
   # Convert to distance format
   d.dist <- as.dist(d.sym)
   # Perform hierarchichal clustering
   h <- hclust(d = d.dist, method = "average")
   # Convert to dendrogram
   dend <- as.dendrogram(h)
   
   
   ####### Calculate groups #######
message("step 3")
   if (!(is.null(ANI_threshold))) {
     clusters <- cutree(tree = dend, h = 100 - ANI_threshold)
     nclus <<- length(table(clusters))
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
     cluster_assignments <<- data.frame(allgenomes.cluster, 
                                       allgenomes.genome, allgenomes.medoid)
     colnames(cluster_assignments) <<- c("Cluster", "Genome", 
                                        "Cluster_Medoid")
     
     clusters_order <- cutree(tree = dend, h = 100 - ANI_threshold, order_clusters_as_data = FALSE)
     clusters_order <- as.data.frame(clusters_order)
   
   ###### Bootstraping dendrogram #####
message("step 4")
   message(paste("Bootstrap resampling", bootstrap_iterations, "times ..."))
   bs.results <<- clusterboot(d.dist, B=bootstrap_iterations,
                       distances=TRUE,
                       k=nclus,
                       multipleboot=TRUE,
                       bootmethod = "boot",
                       clustermethod = disthclustCBI,
                       method="average",
                       count=FALSE,
                       seed = seed)
   
   #print(bs.results$bootmean)
     
   ###### Make dendrogram #######
   dend <<- dend %>% set("labels_cex", label_size)

   # Return clusters
   return(list(medoid_genomes = medoid_genomes,
               cluster_assignments = cluster_assignments))
   }
   else {
     return(NULL)
   }
 }
 
 # Plot dendrograms custom function
 plot.dendrogram.custom <- function(xline = NULL, xlinecol = "#20A387FF", 
                                    xlinetype = "dashed"){
   par(mar = c(10,10,10,10))
   
   
   cluster_ids <- split(as.data.frame(bs.results$partition), bs.results$partition) # split vector of clusters groups
   xy <- get_nodes_xy(dend) # get all coordinates for nodes
   
   
   plot(dend, horiz = T, xlab = "Dissimilarity (100-ANI)", 
        main = taxonomy, xlim=c(25,0), cex.main=0.8) %>%
     abline(v = xline, lty = xlinetype, col = xlinecol)
   
   # accessing the ids of the elements per cluster
   for (name in names(cluster_ids)) {
     cluster_elements <- row.names(cluster_ids[[name]])
     
     # THINK ABOUT IMPLEMENTING HERE
     parent_nodes_xy <- which_node(dend, cluster_elements, max_id=FALSE)
     # closest_node_xy <- tail(parent_nodes_xy, n=2)[1] #trying to align to the middle of the edge
     current_node_xy <- max(parent_nodes_xy)
     
     bsvalue <- round(bs.results$bootmean[as.numeric(name)]*100, digits=0)
     # text((xy[current_node_xy,][2] + xy[closest_node_xy,][2] / 2),
     #      xy[current_node_xy,][1],
     #      labels=bsvalue, 
     #      cex=0.7, font=2,
     #      col=ifelse(bsvalue >= 85,"blue","red"))
     text((xy[current_node_xy,][2]+0.5),# + xy[closest_node_xy,][2] / 2),
          xy[current_node_xy,][1]+1.2,
          labels=bsvalue,
          cex=0.7, font=2,
          col=ifelse(bsvalue >= 85,"blue","red"))
   }
   
   #print(cluster_assignments$Genome)
   
   dend %>% color_branches(dend, k=nclus, groupLabels = T) %>% color_labels(dend, k=nclus) %>%
   #dend %>% color_branches(dend, clusters = cluster_assignments$Cluster, groupLabels = T) %>% color_labels(dend, k=nclus) %>%
     plot(horiz = T, xlab = "Dissimilarity (100-ANI)", main = taxonomy, xlim=c(25,0), cex.main=0.8) %>%
     abline(v = xline, lty = xlinetype, col = xlinecol)
 }
 
###### SETTING DIRECTORY AND OTHER FILES #####

# Define whole path
setwd(directory)

# Define current directory
current_dir <- tail(unlist(strsplit(getwd(), "/")), n=1)

# Generate plots file
pdf(paste0("plots_",current_dir,".pdf"), width = 10, height = 10)

# Read GTDB taxonomy file
taxonomy <- readChar(taxonomy_file, file.info(taxonomy_file)$size)

# Read fastANI output file
ani <- read.ANI.custom(file = fastani_output_file)

comparisons <- length(ani$ANI)
total_genomes <- length(unique(ani$query))

###### DEFINING OTHER FUNCTIONS ######

# calculate dendrogram function
calculate_dendrogram <- function(){
  # Generate plot and return clusters
  cluster_dendrogram <- ANI.dendrogram.custom(ani_object = ani,
                                              ANI_threshold = ani_threshold,
                                              label_size = 0.5)
  
  # Write number of genome per clusters to file
  counts <- as.data.frame(table(cluster_dendrogram$cluster_assignments$Cluster))
  colnames(counts) <- c("Cluster","Counts")
  write.table(counts,
              file='cluster_counts.tsv',
              quote=FALSE, sep='\t', 
              row.names=FALSE)
  
  # Write clusters to file with bootstrap values
  clusters_plus_bs <- cbind(cluster_dendrogram$cluster_assignments,
                            round(bs.results$bootmean[cluster_dendrogram$cluster_assignments$Cluster]*100, digits=0))
  names(clusters_plus_bs)[4] <- "BootstrapValue"
  
  #write.table(cluster_dendrogram$cluster_assignments), 
  write.table(clusters_plus_bs,
              file='cluster_summary.tsv', 
              quote=FALSE, sep='\t', 
              row.names=FALSE)
  
  # Write total number of clusters to file
  cluster_number <<- length(cluster_dendrogram$medoid_genomes$Cluster)
  write(cluster_number, file='total_clusters.txt')
  
}

# plot paiwise distances histogram function
plot_histogram <- function(){
  ggplot(ani, aes(ANI)) +
  geom_histogram(binwidth = 0.01) +
    scale_x_continuous(breaks=seq(75,100,1)) +
    ggtitle(paste(comparisons,"pairwise comparisons", "(", total_genomes, "x", total_genomes, ")")) + xlab("ANI distance") + ylab("Count") +
    labs(subtitle = paste("GTDB Taxonomy:",taxonomy),
         caption = "Note: fastANI is made to work with ANI in the range of 80-100")
}
  
# plot final summary
plot_final_summary <- function(){
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.6, paste("Summary\n",current_dir), cex = 2.5, col = "black")
  text(x = 0.5, y = 0.4, paste("GTDB Taxonomy:\n",taxonomy), cex = 0.8, col = "black")
  text(x = 0.5, y = 0.5, paste("Number of bins:\n",total_genomes), cex = 1.0, col = "black")
  text(x = 0.5, y = 0.3, paste0("Number of clusters"," (ANI ",ani_threshold, "):", "\n",cluster_number), cex = 2.0, col = "red")
}

###### CALLING FUNCTIONS ######
message("Calculating clusters ...")
message("Writting summary files ...")
calculate_dendrogram()

message("Plotting summary ...")
plot_final_summary()

message("Plotting dendrogram ...")
plot.dendrogram.custom(xline = 100-ani_threshold,
                       xlinecol = "#20A387FF",
                       xlinetype = "dashed")

message("Plotting histogram ...")
plot_histogram() # pairwise comparisons histogram

# Plot blank canvas
# plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
