#!/usr/bin/env Rscript
#Matti Ruuskanen, Nov 2016
#Combines swarm output clusters to matching OTU table to sum the centroid totals

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE, lib = "~/.local/lib64/R/library/")
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("optparse", "parallel", "doParallel", "foreach")
ipak(packages)

ncores <- detectCores()

option_list <- list(
  make_option(
    c("-i", "--swarm_file"),
    type = "character",
    default = NULL,
    help = "swarm output file with clusters on rows",
    metavar = "character"
  ),
  make_option(
    c("-t", "--otu_table"),
    type = "character",
    default = NULL,
    help = "otu table of dereplicated sequences",
    metavar = "character"
  ),
  make_option(
    c("-p", "--processors"),
    type = "numeric",
    default = ncores,
    help = "number of cores to be used [default = all cores]",
    metavar = "numeric"
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    default = getwd(),
    help = "output folder [default = %default]",
    metavar = "character"
  )
)


opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)

if (is.null(opt$swarm_file) | is.null(opt$otu_table)) {
  print_help(opt_parser)
  stop("Both swarm output table and otu table with dereplicated sequences must be supplied.", call. = FALSE)
}

setwd(opt$out)

swarm <- readLines(opt$swarm_file)
swarm <- strsplit(swarm, "[[:space:]]+")
names(swarm) <- sapply(swarm, `[[`, 1)
names(swarm) <- sub("(.*?_[0-9]*?)_.*", "\\1", names(swarm))
swarm <- lapply(swarm, function(x) sub("(.*?_[0-9]*?)_.*", "\\1", x))

otu_table <- read.table(opt$otu_table, header=T)

centroids <- list()
list <- list()
dirs <- "centroids"
for (i in dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
  }
}

targetcentroids <- names(swarm)
targetcentfiles <- paste(targetcentroids, ".csv", sep = "")
todocent <- targetcentroids[!targetcentfiles %in% dir(dirs[1])]
parresult <-
  foreach (i = unique(todocent), .combine = rbind) %dopar% {
    otus <- unlist(swarm[names(swarm) %in% i])
    out <- as.data.frame(otu_table[which(otu_table$OTU %in% otus),])
    out <- t(data.frame(colSums(out[,-1])))
    out <- cbind(OTU = i, out)
      write.csv(out, paste("centroids/", i, ".csv", sep = ""), row.names=F)
    }

j <- 0
for (i in dir("centroids")) {
  j <- j + 1
  centroids[[j]] <- read.csv(paste("centroids/", i, sep = ""))
}

df <- do.call(rbind.data.frame, centroids)

write.table(
  df,
  file = paste(opt$out, "/swarm_otu_table.txt", sep = ""),
  row.names = F,
  quote = F
)

unlink("centroids", recursive = T)
