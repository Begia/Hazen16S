# # installing/loading the package:
# if(!require(installr)) {
#   install.packages("installr"); require(installr)} #load / install+load installr
# 
# # using the package:
# updateR() # this will start the updating process of your R installation.  It will check for newer versions, and if one is available, will guide you through the decisions you'd need to make.

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 

packages <- c("biomformat", "ggplot2", "phangorn", "ape", "plyr", "magrittr", "scales", "data.table", "ranger",
              "knitr", "edarf", "checkmate", "hopach", "svglite", "caret", "e1071",
              "ggthemes", "ggfortify", "devtools", "phyloseq", "vegan", "stringi", "lme4", "quantreg", "ggbiplot", "NbClust", "Rtsne", "dbscan")

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, dependencies = T, suppressUpdates = T)
  sapply(pkg, require, character.only = TRUE)
}

ipak(packages)


#read in functions from miseq (https://github.com/michberr/MicrobeMiseq/blob/master/R/miseqR.R)
source("C:/Users/Matti/Dropbox/Matti-Alex/R/miseqR.R")


theme_set(theme_tufte(base_family = "sans", base_size = 18) + theme(panel.border = element_rect(colour = "black", fill = NA), 
                                                    axis.text = element_text(colour = "black", size = 18)))

#function for calculating standard error
se <- function(x) sqrt(var(x)/length(x))

#veganifyOTU function extracted from internal phyloseq function
veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#This is an import of a corrected rarefaction curve function 
#from "http://stackoverflow.com/questions/22714775/how-do-i-colour-lines-separately-in-rarecurve-vegan-package"
rarec <- function (x, step = 20, sample, xlab = "Sample Size", ylab = "OTUs", 
                   label = TRUE, cols = c(rep('red', nrow(x))), ...) {
  tot <- rowSums(x)
  S <- specnumber(x)
  nr <- nrow(x)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_len(length(out))) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = cols[ln], ...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

#import the tree, normalized biom file, mapping files and SINA metadata
hazen_tree <- read_tree("D:/VirtualBox/VirtualBox Share/16S/Hazen16S/hazen_tree.tre")
hazensummer_a_tree <- read_tree("D:/VirtualBox/VirtualBox Share/16S/hazensummer/hazensummer_a_tree.tre")
hazensummer_b_tree <- read_tree("D:/VirtualBox/VirtualBox Share/16S/hazensummer/hazensummer_b_tree.tre")
trees <- list(hazen_tree, hazensummer_a_tree, hazensummer_b_tree)

hazen_all_biom <- import_biom("D:/VirtualBox/VirtualBox Share/16S/Hazen16S/otu_table_normalized.biom")
hazensummer_a_all_biom <- import_biom("D:/VirtualBox/VirtualBox Share/16S/hazensummer/a_otu_table_normalized.biom")
hazensummer_b_all_biom <- import_biom("D:/VirtualBox/VirtualBox Share/16S/hazensummer/b_otu_table_normalized.biom")
sample_names(hazen_all_biom) <- substring(sample_names(hazen_all_biom),2)
sample_names(hazensummer_a_all_biom) <- substring(sample_names(hazensummer_a_all_biom),2)
sample_names(hazensummer_b_all_biom) <- substring(sample_names(hazensummer_b_all_biom),2)
norm_bioms <- list(hazen_all_biom, hazensummer_a_all_biom, hazensummer_b_all_biom)

hazen_all_mapping <- import_qiime_sample_data("D:/VirtualBox/VirtualBox Share/16S/Hazen16S/hazen_map_meta_corrected.txt")
hazensummer_a_all_mapping <- import_qiime_sample_data("D:/VirtualBox/VirtualBox Share/16S/hazensummer/hazensummer_a_map_meta_corrected.txt")
hazensummer_b_all_mapping <- import_qiime_sample_data("D:/VirtualBox/VirtualBox Share/16S/hazensummer/hazensummer_b_map_meta_corrected.txt")
mapping_files <- list(hazen_all_mapping, hazensummer_a_all_mapping, hazensummer_b_all_mapping)

SINA_csv <- read.csv("D:/VirtualBox/VirtualBox Share/16S/Hazen16S/swarm_cluster_seeds_rna_aligned.fasta.csv",stringsAsFactors = F)
hazensummer_a_SINA_csv <- read.csv("D:/VirtualBox/VirtualBox Share/16S/hazensummer/a_swarm_cluster_seeds_rna_aligned.fasta.csv",stringsAsFactors = F)
hazensummer_b_SINA_csv <- read.csv("D:/VirtualBox/VirtualBox Share/16S/hazensummer/b_swarm_cluster_seeds_rna_aligned.fasta.csv",stringsAsFactors = F)
to_refine <- list(SINA_csv, hazensummer_a_SINA_csv, hazensummer_b_SINA_csv)

dataset_names <- c("hazen_all", "hazensummer_a", "hazensummer_b")

#refine taxonomy of hits >99% id to species level
SINA_refined <- list(NULL)
tax_tables <- list(NULL)
phyloseqs <- list(NULL)
SILVA_128_taxonomy <- readLines("D:/VirtualBox/VirtualBox Share/16S/SILVA_128.taxonomy")
SILVA_128_taxonomy <- data.frame(acc=sub("^(.*?) {1}.*","\\1",SILVA_128_taxonomy),
                                 taxonomy=sub("^.*? {1}(.*?)","\\1",SILVA_128_taxonomy))
omnipresent_otus <- list(NULL)
for (i in 1:length(to_refine)) {
over_99_rownumbers <- as.numeric(rownames(to_refine[[i]][to_refine[[i]]$align_ident_slv > 99,]))
over_99_id <- to_refine[[i]][over_99_rownumbers,15]
over_99_id <- sub("^(.*?)~.*","\\1",over_99_id)
over_99_id <- sub("\\.[0-9]","",over_99_id)
over_99 <- data.frame(acc=over_99_id,rownumber=over_99_rownumbers)
tax_refinement <- SILVA_128_taxonomy[match(over_99$acc, SILVA_128_taxonomy$acc, nomatch=0),]
over_99 <- merge(tax_refinement,over_99,by="acc")
SINA_refined[[i]] <- to_refine[[i]]
SINA_refined[[i]][match(over_99$rownumber,rownames(SINA_refined[[i]])),] <- over_99$taxonomy
SINA_refined[[i]] <- SINA_refined[[i]][,14]

#create tax tables
tax <- strsplit(SINA_refined[[i]],";")
tax <- as.matrix(rbind.fill(lapply(tax,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)})))
rownames(tax) <- to_refine[[i]]$name
#if an 8th tax level exists (as Eukaryotes sometimes have), remove it
if (ncol(tax) > 7) {tax <- tax[,1:7]}
colnames(tax) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_tables[[i]] <- tax_table(tax)

#construct phyloseq objects
phyloseqs[[i]] <- phyloseq(norm_bioms[[i]], tax_tables[[i]], mapping_files[[i]])

#remove unclassified sequences, chloroplasts and mitochondria (and bacteria from archaeal data set and archaea from bacterial data set)
phyloseqs[[i]] <- subset_taxa(phyloseqs[[i]], !(Domain %in% "Unclassified"))
phyloseqs[[i]] <- subset_taxa(phyloseqs[[i]], !(Domain %in% "Eukaryota"))
phyloseqs[[i]] <- subset_taxa(phyloseqs[[i]], !(Class %in% "Chloroplast"))  
phyloseqs[[i]] <- subset_taxa(phyloseqs[[i]], !(Family %in% "Mitochondria"))
if (i == 2) {
  phyloseqs[[i]] <- subset_taxa(phyloseqs[[i]], Domain %in% "Archaea")
}
if (i == 3) {
  phyloseqs[[i]] <- subset_taxa(phyloseqs[[i]], Domain %in% "Bacteria")
}

phyloseqs[[i]] <- phyloseq(otu_table(phyloseqs[[i]]), tax_tables[[i]], mapping_files[[i]], trees[[i]])
}
#garbage collection, free up RAM
rm(SILVA_128_taxonomy, tax_refinement, over_99, to_refine, SINA_refined, SINA_csv, hazensummer_a_SINA_csv, hazensummer_b_SINA_csv)
gc()


#set sample levels for the main and all following data sets
sample_data(phyloseqs[[1]])$X.SampleID <- paste(sub(".*([0-9]{4})$","\\1",sample_data(phyloseqs[[1]])$Description), gsub("\\.", " ", sample_data(phyloseqs[[1]])$Site),sample_data(phyloseqs[[1]])$Depth.cm,"cm",sep=" ",collapse=NULL)
sample_data(phyloseqs[[2]])$X.SampleID <- paste(gsub("\\.", " ", sample_data(phyloseqs[[2]])$Site),sample_data(phyloseqs[[2]])$Depth.cm,"cm",sep=" ",collapse=NULL)
sample_data(phyloseqs[[3]])$X.SampleID <- paste(gsub("\\.", " ", sample_data(phyloseqs[[3]])$Site),sample_data(phyloseqs[[3]])$Depth.cm,"cm",sep=" ",collapse=NULL)

#remove John's Island
#phyloseqs[[1]] <- subset_samples(phyloseqs[[1]], !(Site %in% "Johns.Island"))
#phyloseqs[[1]] <- prune_taxa(taxa_sums(phyloseqs[[1]]) > 0, phyloseqs[[1]])

sample_levels <-
  c("2014 Snowgoose Bay 0.25 cm",
    "2014 Snowgoose Bay 0.5 cm",
    "2014 Snowgoose Bay 0.75 cm",
    "2014 Snowgoose Bay 1 cm",
    "2014 Snowgoose Bay 5.5 cm",
    "2014 Snowgoose Bay 8 cm",
    "2014 Snowgoose Bay 27.5 cm",
    "2015 Snowgoose Bay 1 cm",
    "2015 Snowgoose Bay 2 cm",
    "2015 Snowgoose Bay 3 cm",
    "2015 Snowgoose Bay 4 cm",
    "2015 Snowgoose Bay 5 cm",
    "2014 Deep Hole 0.25 cm",
    "2014 Deep Hole 0.5 cm",
    "2014 Deep Hole 0.75 cm",
    "2014 Deep Hole 1 cm",
    "2014 Deep Hole 5.5 cm",
    "2014 Deep Hole 8 cm",
    "2015 Deep Hole 1 cm",
    "2015 Deep Hole 2 cm",
    "2015 Deep Hole 3 cm",
    "2015 Deep Hole 4 cm",
    "2015 Deep Hole 5 cm",
    "2015 Johns Island 1 cm",
    "2015 Johns Island 2 cm",
    "2015 Johns Island 3 cm",
    "2015 Johns Island 4 cm",
    "2015 Johns Island 5 cm",
    "2015 Skeleton Lake 1 cm",
    "2015 Skeleton Lake 2 cm",
    "2015 Skeleton Lake 3 cm",
    "2015 Skeleton Lake 4 cm",
    "2015 Skeleton Lake 5 cm"
  )

sample_levels_summer <-
  c("Pond1 0.5 cm",  
    "Pond1 1 cm",  
    "Pond1 1.5 cm", 
    "Pond1 2 cm",  
    "Pond1 2.5 cm",  
    "Pond1 3 cm",  
    "Pond1 3.5 cm",  
    "Pond1 4 cm",
    "Skeleton Lake 0.5 cm", 
    "Skeleton Lake 1 cm",
    "Skeleton Lake 1.5 cm",
    "Skeleton Lake 2 cm", 
    "Skeleton Lake 2.5 cm",
    "Skeleton Lake 3 cm",
    "Skeleton Lake 3.5 cm",
    "Skeleton Lake 4 cm", 
    "Skeleton Lake 4.5 cm" ,
    "Skeleton Lake 5 cm", 
    "Skeleton Lake 5.5 cm", 
    "Skeleton Lake 6 cm"
  )

sample_data(phyloseqs[[1]])$X.SampleID <- factor(sample_data(phyloseqs[[1]])$X.SampleID, levels = sample_levels)
sample_data(phyloseqs[[2]])$X.SampleID <- factor(sample_data(phyloseqs[[2]])$X.SampleID, levels = sample_levels_summer)
sample_data(phyloseqs[[3]])$X.SampleID <- factor(sample_data(phyloseqs[[3]])$X.SampleID, levels = sample_levels_summer)

#remove samples with lacking metadata and unsuited for comparisons (deep sediment layers from 2014)
phyloseqs[[1]] <- subset_samples(phyloseqs[[1]], Depth.cm <= 5.0)
phyloseqs[[1]] <- prune_taxa(taxa_sums(phyloseqs[[1]]) > 0, phyloseqs[[1]])
sample_levels <- c(levels(sample_data(phyloseqs[[1]])$X.SampleID))

for (i in 1:length(phyloseqs)) {
sample_names(phyloseqs[[i]]) <- sample_data(phyloseqs[[i]])$X.SampleID

#at this point the data is ready for analysis

#write quality controlled biom files for FAPROTAX
tax_table <- data.frame(tax_table(phyloseqs[[i]]))
tax0 <- do.call("paste", c(list(rownames(tax_table)), tax_table, sep="; "))
tax0 <- sub("^.+(Bacteria.*?); NA.*", "\\1", tax0)
tax0 <- sub("^.+(Archaea.*?); NA.*", "\\1", tax0)
tax0 <- sub("^OTU_[0-9]*; ", "", tax0)
otu_table <- data.frame(otu_table(phyloseqs[[i]]))
colnames(otu_table) <- sample_names(phyloseqs[[i]])
otu_table$taxonomy <- tax0
otu_table <- cbind(OTU=rownames(otu_table), otu_table)
if (i == 1) {
  write.table(otu_table, "D:/VirtualBox/VirtualBox Share/16S/Hazen16S/hazen_all_qcd.txt", quote=F, sep="\t", row.names = F)
} else {
  write.table(otu_table, paste0("D:/VirtualBox/VirtualBox Share/16S/hazensummer/", dataset_names[i], "_qcd.txt"), quote=F, sep="\t", row.names = F)
}
}

#the following commands ran at CAC cluster
# biom convert -i ~/matti/Hazen16S/hazen_all_qcd.txt -o ~/matti/Hazen16S/hazen_all_qcd.biom --table-type="OTU table" --to-json
# python ~/bin/FAPROTAX_1.0/collapse_table.py -i ~/matti/Hazen16S/hazen_all_qcd.biom -o ~/matti/Hazen16S/hazen_func_table.biom -g ~/bin/FAPROTAX_1.0/FAPROTAX_Hazen.txt --collapse_by_metadata 'taxonomy' --group_leftovers_as 'other' --out_group_overlaps ~/matti/Hazen16S/hazen_func_table_overlaps.txt --output_format_group_overlaps classical  -l ~/matti/Hazen16S/hazen_FAPROTAX.log --disable_group_set_operations --out_groups2records_table ~/matti/Hazen16S/hazen_func_table_groups.txt -v --force
# biom convert -i ~/matti/Hazen16S/hazen_func_table.biom -o ~/matti/Hazen16S/hazen_func_table.txt --to-tsv
# biom convert -i ~/matti/hazensummer_a/hazensummer_a_qcd.txt -o ~/matti/hazensummer_a/hazensummer_a_qcd.biom --table-type="OTU table" --to-json
# python ~/bin/FAPROTAX_1.0/collapse_table.py -i ~/matti/hazensummer_a/hazensummer_a_qcd.biom -o ~/matti/hazensummer_a/hazensummer_a_func_table.biom -g bin/FAPROTAX_1.0/FAPROTAX_Hazen.txt --collapse_by_metadata 'taxonomy' --group_leftovers_as 'other' --out_group_overlaps ~/matti/hazensummer_a/hazensummer_a_func_table_overlaps.txt --output_format_group_overlaps classical -l ~/matti/hazensummer_a/hazensummer_a_FAPROTAX.log --disable_group_set_operations --out_groups2records_table ~/matti/hazensummer_a/hazensummer_a_func_table_groups.txt -v --force
# biom convert -i ~/matti/hazensummer_a/hazensummer_a_func_table.biom -o ~/matti/hazensummer_a/hazensummer_a_func_table.txt --to-tsv
# biom convert -i ~/matti/hazensummer_b/hazensummer_b_qcd.txt -o ~/matti/hazensummer_b/hazensummer_b_qcd.biom --table-type="OTU table" --to-json
# python ~/bin/FAPROTAX_1.0/collapse_table.py -i ~/matti/hazensummer_b/hazensummer_b_qcd.biom -o ~/matti/hazensummer_b/hazensummer_b_func_table.biom -g bin/FAPROTAX_1.0/FAPROTAX_Hazen.txt --collapse_by_metadata 'taxonomy' --group_leftovers_as 'other' --out_group_overlaps ~/matti/hazensummer_b/hazensummer_b_func_table_overlaps.txt --output_format_group_overlaps classical -l ~/matti/hazensummer_b/hazensummer_b_FAPROTAX.log --disable_group_set_operations --out_groups2records_table ~/matti/hazensummer_a/hazensummer_a_func_table_groups.txt -v --force
# biom convert -i ~/matti/hazensummer_b/hazensummer_b_func_table.biom -o ~/matti/hazensummer_b/hazensummer_b_func_table.txt --to-tsv

#further conversion and FAPROTAX analysis done here in one big for loop running once for each separate dataset
bioms <- list(NULL)
modelstorun <- list(NULL)
richness_list <- list(NULL) 
diversity_list <- list(NULL)
phyla_abundances <- list(NULL)
phyla_relative_abundances <- list(NULL)
func_abundances <- list(NULL)
func_relative_abundances <- list(NULL)
mantels <- list(NULL)
func_data <- list(NULL)
for (i in 1:length(phyloseqs)) {
#reimport the FAPROTAX data
if (i == 1) {
  func <- read.csv("D:/VirtualBox/VirtualBox Share/16S/Hazen16S/hazen_func_table.txt",sep="\t",skip=1,row.names=1)
  colnames(func) <- sub("^.", "", colnames(func))
} else {
  func <- read.csv(paste0("D:/VirtualBox/VirtualBox Share/16S/hazensummer/", dataset_names[i], "_func_table.txt"),sep="\t",skip=1,row.names=1)
}

colnames(func) <- sample_names(phyloseqs[[i]])
func <- phyloseq(otu_table(func,taxa_are_rows = T), sample_data(phyloseqs[[i]]))
func <- prune_taxa(taxa_sums(func)>0, func)
func_pruned <- prune_taxa(!(taxa_names(func) %in% "other"),func)



#capitalize function names and remove underscores
rownames(otu_table(func_pruned)) <- paste(toupper(substr(rownames(otu_table(func_pruned)), 1, 1)), substr(rownames(otu_table(func_pruned)), 2, nchar(rownames(otu_table(func_pruned)))), sep="")
rownames(otu_table(func_pruned)) <- gsub("_", " ", rownames(otu_table(func_pruned)))

#analysis producing graphs/results starts from here
#make rarefaction curves:
#read in non-normalized biom table
if (i == 1) {
  biom_non_normalized <- import_biom("D:/VirtualBox/VirtualBox Share/16S/Hazen16S/swarm_otu_table_with_singletons.biom")
} else {
  biom_non_normalized <- import_biom(paste0("D:/VirtualBox/VirtualBox Share/16S/hazensummer/", dataset_names[i], "_swarm_otu_table_with_singletons.biom"))
}
#remove samples that weren't analyzed from the 2014/2015 hazen dataset
if (i == 1) {biom_non_normalized <- biom_non_normalized[,-c(5,6,11:13)]}
#correct sample names
sample_names(biom_non_normalized) <- sample_names(phyloseqs[[i]])



#transpose data
data <- t(biom_non_normalized)
#use vegan functions for...

#getting number of species
S <- specnumber(data)
#get the rarefaction cutoff value (minimum number of reads per sample)
raremax <- min(rowSums(data))
#rarefy at the cutoff
Srare <- rarefy(data, raremax)
#plot the effect of rarefaction (rarefied number of species as a function of total observed number of species)
plot(S, Srare, xlab = "Observed No. of OTUs", ylab = "Rarefied No. of OTUs")
abline(0, 1)
#plot rarefaction curves with the imported rarefaction curves "rarec" function
svg(filename=paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_rarefaction.svg"), width=11.8, height=7.8, pointsize=12)
rarec(data, sample = raremax)
dev.off()


#plot taxonomic distribution at a phylum level

#agglomerate taxa to phylum levels and set unknown phyla (NA) to "Unknown"
unknown_OTU <- rownames(tax_table(phyloseqs[[i]])[tax_table(phyloseqs[[i]])[,2] %in% NA,2])
phyla <- merge_taxa(phyloseqs[[i]], unknown_OTU)
#subdivide Proteobacteria to Class level (for the taxa which have a class level)
if(any(tax_table(phyla)[,2] %in% "Proteobacteria")) {
  replacement_proteobacteria <- tax_table(phyla)[tax_table(phyla)[,2] %in% "Proteobacteria",3]
  replacement_proteobacteria[replacement_proteobacteria %in% NA] <- "Proteobacteria"
  tax_table(phyla)[tax_table(phyla)[,2] %in% "Proteobacteria",2] <- replacement_proteobacteria
  }
#merge phyla that have less than 1% abundance (or unknown Phylum) in the data set to "other / unknown"
phyla <- tax_glom(phyla, taxrank = "Phylum")
least_abundant_phyla <- names(taxa_sums(phyla))[taxa_sums(phyla)/sum(taxa_sums(phyla))*100 < 1]
phyla <- merge_taxa(phyla, least_abundant_phyla, 1)
tax_table(phyla)[tax_table(phyla)[,2] %in% NA,2] <- "Other / Unknown"
#plot % abundance of the most abundant Phyla
phyla <- transform_sample_counts(phyla, function(x) 100 * x/sum(x))

#set order of phyla
phyla_levels_all <-
  c("AC1",
    "Acetothermia",
    "Acidobacteria",
    "Actinobacteria",
    "Aminicenantes",
    "Armatimonadetes",
    "Atribacteria",
    "Bacteroidetes",
    "Bathyarchaeota",
    "Candidatus Berkelbacteria",
    "BJ-169",
    "BRC1",
    "Caldiserica",
    "Chlamydiae",
    "Chlorobi",
    "Chloroflexi",
    "CPR2",
    "Cyanobacteria",
    "Deinococcus-Thermus",
    "Elusimicrobia",
    "Euryarchaeota",
    "Miscellaneous Euryarchaeotic Group(MEG)",
    "FBP",
    "FCPU426",
    "Fibrobacteres",
    "Firmicutes",
    "GAL15",
    "Gemmatimonadetes",
    "Gracilibacteria",
    "Hydrogenedentes",
    "Ignavibacteriae",
    "Latescibacteria",
    "LCP-89",
    "Lentisphaerae",
    "Microgenomates",
    "Nitrospinae",
    "Nitrospirae",
    "Omnitrophica",
    "Parcubacteria",
    "PAUC34f",
    "Peregrinibacteria",
    "Planctomycetes",
    "Proteobacteria",
    "Alphaproteobacteria",
    "Betaproteobacteria",
    "Gammaproteobacteria",
    "Deltaproteobacteria",
    "RBG-1 (Zixibacteria)",
    "Saccharibacteria",
    "Spirochaetae",
    "SR1 (Absconditabacteria)",
    "Tectomicrobia",
    "Thaumarchaeota",
    "TM6 (Dependentiae)",
    "Verrucomicrobia",
    "Woesearchaeota (DHVEG-6)",
    "WS2",
    "WS6",
    "WWE3",
    "Other / Unknown")

#color palettes from http://tools.medialab.sciences-po.fr/iwanthue/
darkcols <- c(
  "#567a54",
  "#6676bf",
  "#b4d68b",
  "#5c9ddf",
  "#a97a3d",
  "#43d7e0",
  "#a67b4f",
  "#4dacd7",
  "#d8c379",
  "#8694c6",
  "#6d803c",
  "#6d9cb9",
  "#d2c089",
  "#70cad1",
  "#aa9676",
  "#75d7bc",
  "#71a2a3",
  "#81c48f",
  "#3d957a",
  "#b0c8a9",
  "#999999"
)


bars <- data.frame(t(otu_table(phyla)))
darkcols <- darkcols[(length(darkcols)+1-ncol(bars)):length(darkcols)]
colnames(bars) <- c(tax_table(phyla)[,2])
phyla_levels <- phyla_levels_all[phyla_levels_all %in% colnames(bars)]
bars <- bars[,phyla_levels]
bars <- melt(as.matrix(bars), varnames=c("Sample", "Phylum"), value.name="Abundance")
bars$Phylum <- factor(bars$Phylum, levels = rev(phyla_levels))
phyla_relative_abundances[[i]] <- bars
if (i == 1) {
  for (a in 1:5) {
  indices <- c(4, 10, 15, 21, 27)
  empty_levels <- c("a", "b", "c", "d", "e")
  sample_levels <- append(sample_levels, empty_levels[a], after = indices[a])
  }
  bars <- rbind(bars, data.frame(Sample = empty_levels, Phylum = levels(bars$Phylum)[1], Abundance = 0))
  
  p1 <- ggplot(bars, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_levels), breaks = bars$Sample[nchar(as.character(bars$Sample))!=1]) + geom_bar(stat="identity") + 
    scale_fill_manual(values=rev(darkcols), name = "Phylum", guide = guide_legend(reverse = T)) + theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
    xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("Taxonomy")
} else {
  if (!any(sample_levels_summer %in% "a")) {sample_levels_summer <- append(sample_levels_summer, "a", after = 8)}
  bars <- rbind(bars, data.frame(Sample = "a", Phylum = levels(bars$Phylum)[1], Abundance = 0))
  p1 <- ggplot(bars, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_levels_summer), breaks = bars$Sample[nchar(as.character(bars$Sample))!=1]) + geom_bar(stat="identity") + 
    scale_fill_manual(values=rev(darkcols), name = "Phylum", guide = guide_legend(reverse = T)) + theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5))+
    xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("Taxonomy")
}


#plot functions from FAPROTAX
#group everything except the 15 most common functions together (if there are that many)
if (length(taxa_names(func_pruned)) > 15){
least_abundant_func <- names(sort(taxa_sums(func_pruned),decreasing=T))[16:length(taxa_sums(func_pruned))]
func_pruned_bar <- merge_taxa(func_pruned, least_abundant_func, 1)
rownames(otu_table(func_pruned_bar))[rownames(otu_table(func_pruned_bar)) %in% names(sort(taxa_sums(func_pruned),decreasing=T))[16]] <- "Other classified"
} else {
  func_pruned_bar <- func_pruned
}
func_pruned_bar <- transform_sample_counts(func_pruned_bar, function(x) 100 * x/sum(x))
bars_func <- data.frame(t(otu_table(func_pruned_bar)))
colnames(bars_func) <- gsub("\\.", " ", colnames(bars_func))
if (any(colnames(bars_func) %in% "Other classified")) {
  bars_func <- cbind(bars_func[-which(colnames(bars_func) %in% "Other classified")], bars_func[which(colnames(bars_func) %in% "Other classified")])
}
funcolors <- c("#999999",
               "#567a54",
               "#6676bf",
               "#b4d68b",
               "#5c9ddf",
               "#a97a3d",
               "#43d7e0",
               "#a67b4f",
               "#4dacd7",
               "#d8c379",
               "#8694c6",
               "#6d803c",
               "#6d9cb9",
               "#d2c089",
               "#70cad1",
               "#aa9676",
               "#75d7bc",
               "#71a2a3",
               "#81c48f",
               "#3d957a",
               "#b0c8a9")

if (!any(colnames(bars_func) %in% "Other classified")) {
  funcolors <- rev(funcolors[2:(length(bars_func)+1)])
} else {
  funcolors <- rev(funcolors[1:(length(bars_func))])
  }
bars_func <- melt(as.matrix(bars_func), varnames=c("Sample", "Function"), value.name="Abundance")
bars_func <- rbind(bars_func[bars_func$Function %in% "Other classified",], bars_func[!bars_func$Function %in% "Other classified",])
bars_func <- bars_func[order(-1:-nrow(bars_func)),]
bars_func$Function <- factor(bars_func$Function, levels = rev(levels(bars_func$Function)))
func_relative_abundances[[i]] <- bars_func
if (i == 1) {
  bars_func <- rbind(bars_func, data.frame(Sample = empty_levels, Function = levels(bars_func$Function)[1], Abundance = 0))
  p2 <- ggplot(bars_func, aes(x=factor(Sample), y=Abundance, fill=factor(Function))) + scale_x_discrete(limits=rev(sample_levels), breaks=bars$Sample[nchar(as.character(bars$Sample))!=1]) + geom_bar(stat="identity") + 
    scale_fill_manual(values=rev(funcolors), name = "Group", guide=guide_legend(reverse=T)) + theme(axis.text.y = element_blank(), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
    ylab("Abundance (%)") + coord_flip() + ggtitle("Functional mapping")
} else {
  bars_func <- rbind(bars_func, data.frame(Sample = "a", Function = levels(bars_func$Function)[1], Abundance = 0))
  p2 <- ggplot(bars_func, aes(x=factor(Sample), y=Abundance, fill=factor(Function))) + scale_x_discrete(limits=rev(sample_levels_summer), breaks=bars$Sample[nchar(as.character(bars$Sample))!=1]) + geom_bar(stat="identity") + 
    scale_fill_manual(values=rev(funcolors), name = "Group", guide=guide_legend(reverse=T)) + theme(axis.text.y = element_blank(), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
    ylab("Abundance (%)") + coord_flip() + ggtitle("Functional mapping")
  }

#plot aerobic and anaerobic functions
aerobic_frame <- data.frame(otu_table(func_pruned))
rownames(aerobic_frame) <- taxa_names(func)[-grep("other", taxa_names(func))]
lookup_aerobic <- read.table("D:/VirtualBox/VirtualBox Share/16S/Hazen16S/FAPROTAX_aerobic_table.txt", header = T, row.names = 1, sep = "\t")
aerobic_frame$aerobic <- rownames(aerobic_frame)
aerobic_frame$aerobic  <- mapvalues(aerobic_frame$aerobic, from = rownames(lookup_aerobic), to = as.vector(lookup_aerobic$aerobic))
aerobic_frame <- aggregate(. ~ aerobic, data = aerobic_frame, sum)
rownames(aerobic_frame) <- aerobic_frame$aerobic
aerobic_frame[1] <- list(NULL)
colnames(aerobic_frame) <- colnames(otu_table(func_pruned_bar))
aerobic_frame_bar <- prop.table(as.matrix(aerobic_frame), margin=2)*100
aerobic_frame_bar <- t(aerobic_frame_bar)
aerobic_frame_bar <- melt(as.matrix(aerobic_frame_bar), varnames=c("Sample", "Aerobic"), value.name="Abundance")
levels(aerobic_frame_bar$Aerobic) <- c("yes", "variable", "no")
aerobiccolors <- c("#2e8286", "#999999", "#6f2f2c")

if (i == 1) {
  image <- ggplot(aerobic_frame_bar, aes(x=factor(Sample), y=Abundance, fill=forcats::fct_rev(factor(Aerobic)))) + scale_x_discrete(limits=sample_levels) + geom_bar(stat="identity") + 
    scale_fill_manual(values=aerobiccolors, name = "Aerobic") + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, size = 14), panel.border = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    xlab("Sample") + ylab("Abundance (%)")
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_aerobic_bars.svg"), plot=image, units="mm", width=300, height=200)
} else {
  image <- ggplot(aerobic_frame_bar, aes(x=factor(Sample), y=Abundance, fill=forcats::fct_rev(factor(Aerobic)))) + scale_x_discrete(limits=sample_levels_summer) + geom_bar(stat="identity") + 
    scale_fill_manual(values=aerobiccolors, name = "Aerobic") + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, size = 14), panel.border = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    xlab("Sample") + ylab("Abundance (%)")
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_aerobic_bars.svg"), plot=image, units="mm", width=300, height=200)
}

bioms[[i]] <- biom_non_normalized
non_normalized <- phyloseq(bioms[[i]], tax_tables[[i]], sample_data(phyloseqs[[i]]), trees[[i]])

#fix sample data at this point
if (i == 1) {
  bad_columns <- c("X.SampleID", "BarcodeSequence", "LinkerPrimerSequence", "BarcodeName", "ReversePrimer")
  column_order <- c("Lake", "Site", "Year", "Description", "H2S.uM", "Depth.cm", "Wdepth.m", "pH", "Redox.mV", "O2.mgL", "n.seq")
  } else {
    # Projectname and Water depth are not needed for the Summer 2015 set
  bad_columns <- c("X.SampleID", "BarcodeSequence", "LinkerPrimerSequence", "BarcodeName", "ProjectName", "ReversePrimer", "Wdepth.m")
  column_order <- c("Site", "Description", "Depth.cm", "pH", "O2.mgL", "NO3.mgL", "Cl.mgL", "SO42.mgL", "n.seq")
}
sample_data_replacement <- data.frame(sample_data(phyloseqs[[i]])[, !(colnames(sample_data(phyloseqs[[i]])) %in% bad_columns)])
sample_data_replacement$Depth.cm <- as.numeric(as.character(sample_data_replacement$Depth.cm))
#get number of sequences as number of non-normalized non-singleton reads per sample
sample_data_replacement$n.seq <- as.numeric(sample_sums(prune_taxa(!taxa_sums(non_normalized) == 1,non_normalized)))
#sample_data_replacement$n.seq <- as.numeric(sample_sums(phyloseqs[[i]]))

if (i == 1) {
  sample_data_replacement <- within(sample_data_replacement, {Lake = ifelse(Site %in% "Skeleton.Lake", "Skeleton.Lake", "Lake.Hazen")})
  sample_data_replacement$Lake <- gsub("\\.", " ", sample_data_replacement$Lake)
  sample_data_replacement$Lake <- factor(sample_data_replacement$Lake)
  sample_data_replacement$ProjectName <- as.vector(sample_data_replacement$ProjectName)
  sample_data_replacement$ProjectName <- sub(".*([0-9]{4})$","\\1",sample_data_replacement$Description)
  sample_data_replacement <- plyr::rename(sample_data_replacement, c("ProjectName" = "Year"))
  sample_data_replacement$Year <- as.factor(sample_data_replacement$Year)
  sample_data_replacement$Site <- as.character(sample_data_replacement$Site)
  sample_data_replacement$Site[sample_data_replacement$Site %in% "Johns.Island"] <- "John's.Island"
} 
sample_data_replacement$Site <- gsub("\\.", " ", sample_data_replacement$Site)
sample_data_replacement$Site <- factor(sample_data_replacement$Site)
sample_data_replacement[apply(sample_data_replacement, 2, function(x) grepl("ND", x))] <- NA
sample_data_replacement <- sample_data_replacement[, apply(sample_data_replacement, 2, function(x) !any(is.na(x)))]

#microprobe data
if (i == 1) {
  microprobe_data <- read.csv("D:/VirtualBox/VirtualBox Share/16S/hazen_all_microprobes.csv",sep="\t")
  microprobe_2014 <- microprobe_data[grep("2014", microprobe_data$Core),]
  microprobe_2014 <- transform(microprobe_2014, bin = cut(Depth.mm, breaks = seq(0,10,by=2.5), include.lowest = T))
  microprobe_2015 <- microprobe_data[grep("2015", microprobe_data$Core),]
  microprobe_2015 <- transform(microprobe_2015, bin = cut(Depth.mm, breaks = seq(0,50,by=10), include.lowest = T))
  microprobe_data <- rbind(microprobe_2014, microprobe_2015)

} else {
  microprobe_data <- read.csv("D:/VirtualBox/VirtualBox Share/16S/hazensummer_microprobes.csv",sep="\t")
  microprobe_data <- transform(microprobe_data, bin = cut(Depth.mm, breaks = seq(0,60,by=5), include.lowest = T))
  chemistry_data <- read.csv("D:/VirtualBox/VirtualBox Share/16S/hazensummer_chemistry.csv",sep="\t", row.names = 1)
  sample_rownames <- rownames(sample_data_replacement)[order(sample_data_replacement$Site, sample_data_replacement$Depth.cm)]
  sample_data_replacement <- merge(sample_data_replacement, chemistry_data, by = c("Site", "Depth.cm"))
  sample_data_replacement <- sample_data_replacement[order(sample_data_replacement$Site, sample_data_replacement$Depth.cm),]
  rownames(sample_data_replacement) <- sample_rownames
}

microprobe_data_means <- ddply(microprobe_data, .(Core,bin), numcolwise(median))[-c(1:3)]
microprobe_data_sds <- ddply(microprobe_data, .(Core,bin), numcolwise(sd))[,-c(1:3)]
microprobe_data_sds$Sample <- rownames(sample_data_replacement)
sample_data_replacement <- cbind(sample_data_replacement, microprobe_data_means)
sample_data_replacement <- sample_data_replacement[column_order]


sample_data(phyloseqs[[i]]) <- sample_data_replacement

non_normalized <- phyloseq(bioms[[i]], tax_tables[[i]], sample_data(phyloseqs[[i]]), trees[[i]])

if (i == 1) {
  sample_data_plot <- sample_data_replacement[,-c(1,4,7,11)]
  sample_data_plot$core <- apply(sample_data_plot[,c(2,1) ], 1, paste, collapse = " ")
  sample_data_plot <- sample_data_plot[,-c(1,2)]
  sample_data_plot$Depth.cm <- factor(sample_data_plot$Depth.cm)
  sample_data_plot$Sample <- rownames(sample_data_plot)

#scale variables for the plot
  sample_data_plot$O2.mgL <- sample_data_plot$O2.mgL*(500/15)
  microprobe_data_sds$O2.mgL <- microprobe_data_sds$O2.mgL*(500/15)
  sample_data_plot$H2S.uM <- sample_data_plot$H2S.uM*(500/200)
  microprobe_data_sds$H2S.uM <- microprobe_data_sds$H2S.uM*(500/200)
  sample_data_plot$pH <- sample_data_plot$pH*(500/15)
  microprobe_data_sds$pH <- microprobe_data_sds$pH*(500/15)

  sample_data_plot <- merge(melt(sample_data_plot), melt(microprobe_data_sds), by = c("Sample", "variable"), all = T)
  sample_data_plot$Depth.cm <- factor(sample_data_plot$Depth.cm, levels = rev(levels(sample_data_plot$Depth.cm)))
  sample_data_plot$core <- factor(sample_data_plot$core, levels = c("2014 Snowgoose Bay", "2015 Snowgoose Bay", "2014 Deep Hole", "2015 Deep Hole", "2015 John's Island", "2015 Skeleton Lake"))
  sample_data_plot$sdmin <- replace(sample_data_plot$value.x-sample_data_plot$value.y, sample_data_plot$value.x-sample_data_plot$value.y < 0, 0)
  sample_data_plot$sdmax <- sample_data_plot$value.x+sample_data_plot$value.y

  p3 <- ggplot(sample_data_plot, aes(x=Depth.cm, y=value.x, group=core)) + geom_point(aes(fill = variable, color = variable, shape = variable), size=7, alpha = 0.5) + facet_grid(core~., scales="free_y", space="free_y") + 
  geom_errorbar(aes(ymin=sdmin, ymax=sdmax, color = variable), width = 0.5) + coord_flip() + scale_y_continuous(breaks = seq(0,500,100), limits = c(0,500)) +
  scale_shape_manual(values = c("\u25CF", "\u25B2", "\u25C6", "\u25A0")) + scale_color_manual(values = c("#cc546d", "#6ca65c", "#000000", "#0ac5dc")) + ggtitle("Geochemistry") + #xlab("Depth from sediment surface (cm)") + 
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.spacing = unit(2, "lines"), strip.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_blank(), 
        legend.position="bottom", plot.margin = unit(c(0.8,0.5,6.2,0.5), "lines"), plot.title = element_text(size=28, hjust = 0.5))

  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_taxonomy_multiplot.svg"), plot=p1, units="in", width=10, height=20)
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_chemistry_multiplot.svg"), plot=p3, units="in", width=5, height=19)
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_function_multiplot.svg"), plot=p2, units="in", width=5.8, height=20)
} else if (i == 2) {
  sample_data_plot <- sample_data_replacement[,-c(2,9)]
  sample_data_plot$Depth.cm <- factor(sample_data_plot$Depth.cm)
  sample_data_plot$Sample <- rownames(sample_data_plot)
  
  #scale variables for the plot
  sample_data_plot$SO42.mgL <- sample_data_plot$SO42.mgL/10 #scale for others 0-15, for sulfate 0-150
  
  sample_data_plot <- merge(melt(sample_data_plot), melt(microprobe_data_sds), by = c("Sample", "variable"), all = T)
  sample_data_plot$Depth.cm <- factor(sample_data_plot$Depth.cm, levels = rev(levels(sample_data_plot$Depth.cm)))
  sample_data_plot$sdmin <- replace(sample_data_plot$value.x-sample_data_plot$value.y, sample_data_plot$value.x-sample_data_plot$value.y < 0, 0)
  sample_data_plot$sdmax <- sample_data_plot$value.x+sample_data_plot$value.y
  
  p3 <- ggplot(sample_data_plot, aes(x=Depth.cm, y=value.x, group=Site)) + geom_point(aes(fill = variable, color = variable, shape = variable), size=7, alpha = 0.5) + facet_grid(Site~., scales="free_y", space="free_y") + 
    geom_errorbar(aes(ymin=sdmin, ymax=sdmax, color = variable), width = 0.5) + coord_flip() + scale_y_continuous(breaks = seq(0,12,2), limits = c(0,13)) +
    scale_shape_manual(values = c("\u25B2", "\u25A0", "\u25BC", "\u25C6", "\u25CF")) + scale_color_manual( values = c("#6ca65c", "#0ac5dc", "#f8766d", "#7792db", "#db77ab")) + ggtitle("Geochemistry") + #xlab("Depth from sediment surface (cm)")  +
    theme(axis.text.y = element_text(hjust = 1, size = 22), panel.spacing = unit(2, "lines"), strip.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_blank(), 
          legend.position="bottom", plot.margin = unit(c(0.8,0.5,6.2,0.5), "lines"), plot.title = element_text(size=28, hjust = 0.5))
  
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_taxonomy_multiplot.svg"), plot=p1, units="in", width=8.8, height=10)
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_chemistry_multiplot.svg"), plot=p3, units="in", width=5, height=10.75)
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_function_multiplot.svg"), plot=p2, units="in", width=5.8, height=10.25)
} else {
  p1 <- p1 + guides(fill=guide_legend(nrow=1, reverse = T))
  p2 <- p2 + guides(fill=guide_legend(nrow=1, reverse = T))
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_taxonomy_multiplot.svg"), plot=p1, units="in", width=8.8, height=10)
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_function_multiplot.svg"), plot=p2, units="in", width=5.8, height=10)
}

rm(p1, p2, biom_non_normalized, data)
gc()

#check collinearity and separation of samples on physicochemical variables with PCA
continuous_variables <- colnames(sample_data_replacement[,!(sapply(sample_data_replacement,is.factor))])
continuous_variables <- continuous_variables[!(continuous_variables %in% "n.seq")]
pca_variables <- capture.output(cat(continuous_variables,sep="+"))
pca_from_sample_data <- prcomp(eval(parse(text = paste0("~", pca_variables))),data=sample_data_replacement, scale=T, center=T)
#plot biplots
if (i != 3) {
if (i == 1) {
hazen_all_colours <- c("#910a80", #deep hole
                       "#5199ff", #john's island
                       "#aac856", #skeleton lake
                       "#a67e48") #snowgoose bay
image <- autoplot(pca_from_sample_data, varname.family = "serif", loadings = T, loadings.label = T, loadings.label.size = 6, loadings.colour = "black", loadings.label.family = "serif",
                  loadings.label.repel = T, loadings.label.colour = "wheat4", groups = sample_data_replacement$Site) + 
  geom_point(aes(colour = sample_data_replacement$Site, shape = sample_data_replacement$Lake, size = sample_data_replacement$Depth.cm)) + 
  geom_point(aes(alpha = as.character(sample_data_replacement$Year)), shape = 3, size = 5) +
  scale_colour_manual(values = hazen_all_colours, name = "Site") + 
  scale_shape_manual(values = c(16, 17), name = "Lake") + scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") +
  guides(colour = guide_legend(override.aes = list(size=5)), shape = guide_legend(override.aes = list(size=5))) +
  xlab(paste0("PC1 (", round(100*summary(pca_from_sample_data)$importance[2],1), "% explained var.)")) + ylab(paste0("PC2 (", round(100*summary(pca_from_sample_data)$importance[5],1), "% explained var.)")) +
  scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)")
} else {
  hazensummer_colours <-c("#229373", #pond1
                          "#aac856") #skeleton lake
  image <- autoplot(pca_from_sample_data, varname.family = "serif", loadings = T, loadings.label = T, loadings.label.size = 6, loadings.colour = "black", loadings.label.family = "serif",
                    loadings.label.repel = T, loadings.label.colour = "wheat4", groups = sample_data_replacement$Site) + 
    geom_point(aes(colour = sample_data_replacement$Site, size = sample_data_replacement$Depth.cm)) + 
    scale_color_manual(values = hazensummer_colours, name = "Site") + 
    guides(colour = guide_legend(override.aes = list(size=5))) + 
    xlab(paste0("PC1 (", round(100*summary(pca_from_sample_data)$importance[2],1), "% explained var.)")) + ylab(paste0("PC2 (", round(100*summary(pca_from_sample_data)$importance[5],1), "% explained var.)")) +
    scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)")
  }
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_PCA_biplot.svg"), plot=image, units="mm", width=300, height=300)
}

avgrichness <- list(NULL)
for (j in 1:10) {
set.seed(42*j)
rarefied <- rarefy_even_depth(non_normalized)

richness <- estimate_richness(rarefied)
Goods.cover <- NULL
for (k in 1:nsamples(rarefied)) {
Goods.cover[k] <- 1-(length(which(otu_table(rarefied)[,k] == 1))/sample_sums(rarefied)[k])
}
richness <- cbind(richness, Goods.cover=Goods.cover)
richness$sample <- sample_names(phyloseqs[[i]])
rownames(richness) <- seq(1+nrow(richness)*(j-1),nrow(richness)*j)
richness$set <- dataset_names[i]
avgrichness[[j]] <- richness
}

richness <- do.call("rbind",avgrichness)
richness <- cbind(richness, Site=gsub(".*\\.(.*\\..*)\\.[0-9]", "\\1", richness$sample))
richness_melt <- richness[which(colnames(richness) %in% c("Chao1", "Shannon", "InvSimpson", "Goods.cover", "sample", "Site"))]
richness_melt <- plyr::rename(richness_melt, replace = c("InvSimpson" = "Simpson's dominance", "Goods.cover" = "Good's coverage"))
richness_melt <- melt(richness_melt)


richness_chao <- aggregate(Chao1 ~ sample, richness, mean)
richness_chao <- cbind(richness_chao, aggregate(Chao1 ~ sample, richness, sd)[2])
richness_shannon <- aggregate(Shannon ~ sample, richness, mean)[2]
richness_shannon <- cbind(richness_shannon, aggregate(Shannon ~ sample, richness, sd)[2])
richness_InvSimpson <- aggregate(InvSimpson ~ sample, richness, mean)[2]
richness_InvSimpson <- cbind(richness_InvSimpson, aggregate(InvSimpson ~ sample, richness, sd)[2])
richness_Goodscover <- aggregate(Goods.cover ~ sample, richness, mean)[2]
richness_Goodscover <- cbind(richness_Goodscover, aggregate(Goods.cover ~ sample, richness, sd)[2])
richness2 <- cbind(richness_chao, richness_shannon, richness_InvSimpson, richness_Goodscover)
colnames(richness2) <- c("sample", "Chao1", "Chao1.sd", "Shannon", "Shannon.sd", "InvSimpson", "InvSimpson.sd", "Goods.cover", "Goods.cover.sd")

diversity_list[[i]] <- richness2

if (i == 1) {
  image <- ggplot(richness_melt, aes(x=factor(sample), y=value)) + scale_x_discrete(limits=sample_levels) + geom_boxplot() + 
    facet_wrap(~variable, scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, size = 8)) + xlab("Sample") + ylab("Value")
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_richness.svg"), plot=image, units="mm", width=400, height=200)
} else {
  image <- ggplot(richness_melt, aes(x=factor(sample), y=value)) + scale_x_discrete(limits=sample_levels_summer) + geom_boxplot() + 
    facet_wrap(~variable, scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, size = 8)) + xlab("Sample") + ylab("Value")
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_richness.svg"), plot=image, units="mm", width=400, height=200)
}


######
predictors_richness <- data.frame(sample_data(phyloseqs[[i]]))
#remove n.seq column
predictors_richness <- predictors_richness[-which(colnames(predictors_richness) %in% "n.seq")]
predictors_richness <- predictors_richness[ , colSums(is.na(predictors_richness)) == 0]

row.names(predictors_richness) <- sample_names(non_normalized)
Shannon <- richness2$Shannon
rf_data_richness <- cbind.data.frame(Shannon, predictors_richness, stringsAsFactors = FALSE)
rf_data_richness$Description <- NULL

#run random forest model
set.seed(42)
classify_richness <- ranger(Shannon ~ ., data = rf_data_richness, num.trees=5000, importance="impurity")

#calculate MSPE estimate with CI for the full model
pred1 <- classify_richness$predictions
y <- rf_data_richness$Shannon
n <- length(y)

# psi is the mean squared prediction error (MSPE) estimate
# sigma2 is the estimate of the variance of the MSPE
psi1 <- mean((y - pred1)^2)
sigma21 <- 1/n * var((y - pred1)^2) 
# 95% CI:
full_MSPE_richness <- c(psi1 - 1.96 * sqrt(sigma21), psi1, psi1 + 1.96 * sqrt(sigma21))

#sort the data by feature importance
importances_richness <- sort(importance(classify_richness), decreasing=T)

#add one feature at a time to the model in order of importance and compare mean MSPEs over the training sets

#reorder the training sets
rf_data_richness_2 <- rf_data_richness[c("Shannon", names(importances_richness))]

comp_models <- list(NULL)
MSPE_richnesses <- list(NULL)
for (k in 1:length(importances_richness)) {
  new_data <- data.frame(Shannon=rf_data_richness_2$Shannon, rf_data_richness_2[seq(2,k+1)])
  set.seed(42)
  comp_models[[k]] <- ranger(Shannon ~ ., data = new_data, num.trees=5000, importance="impurity")
  pred2 <- comp_models[[k]]$predictions
  y2 <- rf_data_richness_2$Shannon
  n2 <- length(y2)
  psi2 <- mean((y2 - pred2)^2)
  sigma22 <- 1/n * var((y2 - pred2)^2) 
  # 95% CI:
  MSPE_richnesses[[k]] <- c(psi2 - 1.96 * sqrt(sigma22), psi2, psi2 + 1.96 * sqrt(sigma22))
}

#find the minimum MSPE and get the feature names
MSPE_richnesses <- do.call("rbind", MSPE_richnesses)
MSPEs <- MSPE_richnesses[,2]
nfeatures <- min(which(MSPEs %in% min(MSPEs)))
feature_names <- names(importances_richness)[1:nfeatures]

#gather, plot and save continuous variable partial dependency plot
pd_richness_2 <- partial_dependence(comp_models[[nfeatures]], vars=feature_names, data=rf_data_richness_2[1:(nfeatures+1)], n=c(25,nrow(rf_data_richness_2)))
pd_continuous_2 <- pd_richness_2[,!(colnames(pd_richness_2) %in% names(Filter(is.factor,pd_richness_2)))]
pd_continuous_2 <- pd_continuous_2[!(apply(pd_continuous_2[-ncol(pd_continuous_2)], 1, function(x) all(is.na(x)))),]
types_cont_2 <- c(rep(colnames(pd_continuous_2)[-ncol(pd_continuous_2)],as.vector(colSums(!is.na(pd_continuous_2))[1:(ncol(pd_continuous_2)-1)])))
pd_continuous_2 <- cbind(pd_continuous_2[ncol(pd_continuous_2)], value = na.omit(unlist(pd_continuous_2[-ncol(pd_continuous_2)])), type=types_cont_2)


image <- ggplot(data = pd_continuous_2, aes(value, Shannon)) + geom_line(size = 1) +
  labs(y="Predicted Shannon Diversity", x="")+ facet_wrap(~type, scales = "free_x") +
  theme(strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_richness_continuous_pd_shannon.svg"), plot=image, units="mm", width=400, height=300)

#gather, plot and save categorical variable partial dependency plot
if (ncol(Filter(is.factor,pd_richness_2)) > 0) {
  
  pd_categorical_2 <- cbind(pd_richness_2$Shannon, Filter(is.factor,pd_richness_2))
  colnames(pd_categorical_2)[1] <- colnames(pd_richness_2)[ncol(pd_richness_2)]
  pd_categorical_2 <- pd_categorical_2[!(apply(data.frame(pd_categorical_2[,-1]), 1, function(x) all(is.na(x)))),]
  types_cat_2 <- c(rep(colnames(pd_categorical_2)[-1],as.vector(colSums(!is.na(pd_categorical_2))[2:ncol(pd_categorical_2)])))
  pd_categorical_2 <- cbind(pd_categorical_2[1], value = na.omit(unlist(pd_categorical_2[-1])), type=types_cat_2)
  
  image <- ggplot(data = pd_categorical_2, aes(x=value, y=Shannon, fill=value)) + geom_bar(stat="identity") +
    labs(y="Predicted Shannon diversity", x="") + facet_wrap(~type, scales = "free_x") + 
    scale_fill_manual(values = c("797368", hazen_all_colours[c(3,1,2,4)])) +
    theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size = 12), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_richness_categorical_pd_shannon.svg"), plot=image, units="mm", width=400, height=300)
}

#gather r2 and other data to a table
richness_frame <- data.frame(dataset=dataset_names[i],
                             type="Shannon",
                             n.before = classify_richness$num.independent.variables, 
                             MSPE.before = round(full_MSPE_richness[2],2),
                             CI.95.before = paste0(round(full_MSPE_richness[1],2),"-",round(full_MSPE_richness[3],2)),
                             r2.before = classify_richness$r.squared,
                             n.after = comp_models[[nfeatures]]$num.independent.variables,
                             MSPE.after = MSPE_richnesses[nfeatures,2],
                             CI.95.after = paste0(round(MSPE_richnesses[nfeatures,1],2),"-",round(MSPE_richnesses[nfeatures,3],2)),
                             r2.after = comp_models[[nfeatures]]$r.squared)
richness_list[[i*10]] <- richness_frame

#Simpson's dominance here
Inv.Simpson <- richness2$InvSimpson
rf_data_richness <- cbind.data.frame(Inv.Simpson, predictors_richness, stringsAsFactors = FALSE)
rf_data_richness$Description <- NULL

#run random forest model
set.seed(42)
classify_richness <- ranger(Inv.Simpson ~ ., data = rf_data_richness, num.trees=5000, importance="impurity")

#calculate MSPE estimate with CI for the full model
pred1 <- classify_richness$predictions
y <- rf_data_richness$Inv.Simpson
n <- length(y)

# psi is the mean squared prediction error (MSPE) estimate
# sigma2 is the estimate of the variance of the MSPE
psi1 <- mean((y - pred1)^2)
sigma21 <- 1/n * var((y - pred1)^2) 
# 95% CI:
full_MSPE_richness <- c(psi1 - 1.96 * sqrt(sigma21), psi1, psi1 + 1.96 * sqrt(sigma21))

#sort the data by feature importance
importances_richness <- sort(importance(classify_richness), decreasing=T)

#add one feature at a time to the model in order of importance and compare mean MSPEs over the training sets

#reorder the training sets
rf_data_richness_2 <- rf_data_richness[c("Inv.Simpson", names(importances_richness))]

comp_models <- list(NULL)
MSPE_richnesses <- list(NULL)
for (k in 1:length(importances_richness)) {
  new_data <- data.frame(Inv.Simpson=rf_data_richness_2$Inv.Simpson, rf_data_richness_2[seq(2,k+1)])
  set.seed(42)
  comp_models[[k]] <- ranger(Inv.Simpson ~ ., data = new_data, num.trees=5000, importance="impurity")
  pred2 <- comp_models[[k]]$predictions
  y2 <- rf_data_richness_2$Inv.Simpson
  n2 <- length(y2)
  psi2 <- mean((y2 - pred2)^2)
  sigma22 <- 1/n * var((y2 - pred2)^2) 
  # 95% CI:
  MSPE_richnesses[[k]] <- c(psi2 - 1.96 * sqrt(sigma22), psi2, psi2 + 1.96 * sqrt(sigma22))
}

#find the minimum MSPE and get the feature names
MSPE_richnesses <- do.call("rbind", MSPE_richnesses)
MSPEs <- MSPE_richnesses[,2]
nfeatures <- min(which(MSPEs %in% min(MSPEs)))
feature_names <- names(importances_richness)[1:nfeatures]

#gather, plot and save continuous variable partial dependency plot
pd_richness_2 <- partial_dependence(comp_models[[nfeatures]], vars=feature_names, data=rf_data_richness_2[1:(nfeatures+1)], n=c(25,nrow(rf_data_richness_2)))
pd_continuous_2 <- pd_richness_2[,!(colnames(pd_richness_2) %in% names(Filter(is.factor,pd_richness_2)))]
pd_continuous_2 <- pd_continuous_2[!(apply(pd_continuous_2[-ncol(pd_continuous_2)], 1, function(x) all(is.na(x)))),]
types_cont_2 <- c(rep(colnames(pd_continuous_2)[-ncol(pd_continuous_2)],as.vector(colSums(!is.na(pd_continuous_2))[1:(ncol(pd_continuous_2)-1)])))
pd_continuous_2 <- cbind(pd_continuous_2[ncol(pd_continuous_2)], value = na.omit(unlist(pd_continuous_2[-ncol(pd_continuous_2)])), type=types_cont_2)

image <- ggplot(data = pd_continuous_2, aes(value, Inv.Simpson)) + geom_line(size= 1) +
  labs(y="Predicted Simpson's dominance", x="") + facet_wrap(~type, scales = "free_x") +
  theme(strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_richness_continuous_pd.svg"), plot=image, units="mm", width=400, height=300)

#gather, plot and save categorical variable partial dependency plot
if (ncol(Filter(is.factor,pd_richness_2)) > 0) {

pd_categorical_2 <- cbind(pd_richness_2$Inv.Simpson, Filter(is.factor,pd_richness_2))
colnames(pd_categorical_2)[1] <- colnames(pd_richness_2)[ncol(pd_richness_2)]
pd_categorical_2 <- pd_categorical_2[!(apply(data.frame(pd_categorical_2[,-1]), 1, function(x) all(is.na(x)))),]
types_cat_2 <- c(rep(colnames(pd_categorical_2)[-1],as.vector(colSums(!is.na(pd_categorical_2))[2:ncol(pd_categorical_2)])))
pd_categorical_2 <- cbind(pd_categorical_2[1], value = na.omit(unlist(pd_categorical_2[-1])), type=types_cat_2)

image <- ggplot(data = pd_categorical_2, aes(x = value, y = Inv.Simpson, fill = value)) + geom_bar(stat="identity") +
  labs(y="Predicted Simpson's dominance", x="") + facet_wrap(~type, scales = "free_x") + 
  scale_fill_manual(values = c(hazen_all_colours, "797368", "#3f268a", "#ffaa00")) +
  theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_richness_categorical_pd.svg"), plot=image, units="mm", width=400, height=300)
}

#gather r2 and other data to a table
richness_frame <- data.frame(dataset=dataset_names[i], 
                             type="Inv.Simpson",
                             n.before = classify_richness$num.independent.variables, 
                             MSPE.before = round(full_MSPE_richness[2],2),
                             CI.95.before = paste0(round(full_MSPE_richness[1],2),"-",round(full_MSPE_richness[3],2)),
                             r2.before = classify_richness$r.squared,
                             n.after = comp_models[[nfeatures]]$num.independent.variables,
                             MSPE.after = MSPE_richnesses[nfeatures,2],
                             CI.95.after = paste0(round(MSPE_richnesses[nfeatures,1],2),"-",round(MSPE_richnesses[nfeatures,3],2)),
                             r2.after = comp_models[[nfeatures]]$r.squared)
richness_list[[i]] <- richness_frame

rm(comp_models)
gc()

#plot sample sequencing depth distribution
image <- ggplot(data.frame(sum = sample_sums(phyloseqs[[i]])), aes(sum)) + 
  geom_histogram(color = "black", fill = "indianred") +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  ylab("Sample count")
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_sequencing_depth_distribution.svg"), plot=image, units="mm", width=300, height=300)

#run envfit to look at the effect of categorical variables (and continuous, but those will be run through ordisurf/gam later on)
sample_data(func_pruned) <- sample_data(phyloseqs[[i]])
sample_data(func_pruned) <- sample_data(func_pruned)[,-which(colnames(sample_data(func_pruned)) %in% "n.seq")]
func_distance <- vegdist(wisconsin(sqrt(veganifyOTU(func_pruned))), distance = "bray")
ordu_func <- ordinate(func_pruned, "NMDS", distance = func_distance)

func_sample_data <- data.frame(sample_data(func_pruned)[,which(colnames(sample_data(func_pruned)) %in% continuous_variables)])
func_env_distance <- vegdist(func_sample_data, method="euclidean")
func_mantel <- mantel(func_distance, func_env_distance, method="pearson", permutations=10000)
ef_func <- envfit(ordu_func,sample_data(func_pruned),permu=10000)
sample_data(func_pruned) <- sample_data(phyloseqs[[i]])
func_data[[i]] <- func_pruned

#prune the datasets from rare taxa (less than 0.1 promille) for subsequent analysis with DPCoA
taxa_sum <- sum(taxa_sums(phyloseqs[[i]]))
taxa_pr <- sort((taxa_sums(phyloseqs[[i]])/taxa_sum)*100,TRUE)
common_names <- names(taxa_pr[taxa_pr>=0.01])
common <- prune_taxa(common_names, phyloseqs[[i]])

#tax_table(phyloseqs[[i]])[rownames(tax_table(phyloseqs[[i]])) %in% omnipresent_otus[[i]],]
sample_data(common) <- sample_data(common)[,-which(colnames(sample_data(common)) %in% "n.seq")]

#midpoint root the tree before running DPCoA
phy_tree(common) = midpoint(phy_tree(common))
common_distance <- DPCoA(common)
ordu_common <- ordinate(common, "NMDS", distance = common_distance$RaoDis)
common_sample_data <- data.frame(sample_data(common)[,which(colnames(sample_data(common)) %in% continuous_variables)])
common_env_distance <- vegdist(common_sample_data, method="euclidean")
common_mantel <- mantel(common_distance$RaoDis, common_env_distance, method="pearson", permutations=10000)
ef_common <- envfit(ordu_common,sample_data(common),permu=10000)
sample_data(common) <- sample_data(phyloseqs[[i]])

# #sanity check with removing skeleton lake samples
# common_hazen_only <- subset_samples(common, Lake %in% "Lake Hazen")
# common_hazen_only <- prune_taxa(taxa_sums(common_hazen_only) > 0, common_hazen_only)
# #midpoint root the tree before running DPCoA
# phy_tree(common_hazen_only) = midpoint(phy_tree(common_hazen_only))
# common_hazen_only_diss <- DPCoA(common_hazen_only)
# ordu_common_hazen_only <- ordinate(common_hazen_only, "NMDS", distance = common_hazen_only_diss$RaoDis)
# ef_common_hazen_only <- envfit(ordu_common_hazen_only,sample_data(common_hazen_only),permu=10000)

# #display a dendrogram of the samples based on DPCoA
# common_hclust <- hclust(common_distance$RaoDis, method="average")
# common_tip_labels <- as(get_variable(common, "Site"), "character")
# plot(as.phylo(common_hclust), show.tip.label = TRUE, tip.color = "black")
# tiplabels(common_tip_labels, col = cols, frame = "none", adj = -0.05, 
#           cex = 0.7)
# 
# func_hclust <- hclust(func_distance, method="average")
# func_tip_labels <- as(get_variable(func_pruned, "Site"), "character")
# plot(as.phylo(func_hclust), show.tip.label = TRUE, tip.color = "black")
# tiplabels(func_tip_labels, col = cols, frame = "none", adj = -0.05, 
#           cex = 0.7)


#clustering analysis on spring 2014/2015 samples
if (i == 1) {
cluster_palette <- c("#03ff00",
                     "#d660cd",
                     "#00e1ff",
                     "#ffb100",
                     "#e65141",
                     "#9800ff")


set.seed(42)
rtsne_common <- data.frame(Rtsne(common_distance$RaoDis, is_distance = T, perplexity = 5)$Y)
rownames(rtsne_common) <- sample_names(common)
rtsne_common_clusters <- hdbscan(rtsne_common, xdist = vegdist(rtsne_common, "euclidean"), minPts = 3)
rtsne_common <- cbind(rtsne_common, cluster = rtsne_common_clusters$cluster, site = sample_data(common)$Site, year = sample_data(common)$Year, depth = sample_data(common)$Depth.cm)
image <- ggplot(data = rtsne_common, aes(x = X1, y = X2)) + stat_ellipse(level = 0.68, aes(color = factor(cluster)), linetype = 2, size = 1.5) +
   geom_point(aes(color = factor(cluster), fill = site, size = depth), shape = 21, stroke = 2) +
   geom_point(aes(alpha = as.character(year)), shape = 3, size=5) +
   scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") + scale_colour_manual(name = "Cluster", values = cluster_palette) +
   scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
                                                                                                       "#5199ff", #john's island
                                                                                                       "#aac856", #skeleton lake
                                                                                                       "#a67e48" #snowgoose bay
   ), name = "Site") +
   guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) + theme(axis.title = element_blank())
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_common_cluster_rtsne2.svg"), plot=image,
      units="mm", width=300, height=300)


#all functional mapping groups, Bray-Curtis dissimilarity
set.seed(42)
rtsne_func <- data.frame(Rtsne(func_distance, is_distance = T, perplexity = 5)$Y)
rtsne_func_clusters <- hdbscan(rtsne_func, xdist = vegdist(rtsne_func, "euclidean"), minPts = 3)
rtsne_func <- cbind(rtsne_func, cluster = rtsne_func_clusters$cluster, site = sample_data(func_pruned)$Site, year = sample_data(func_pruned)$Year, depth = sample_data(func_pruned)$Depth.cm)
image <- ggplot(data = rtsne_func, aes(x = X1, y = X2)) + stat_ellipse(level = 0.68, aes(color = factor(cluster)), linetype = 2, size = 1.5) +
  geom_point(aes(color = factor(cluster), fill = site, size = depth), shape = 21, stroke = 2) +
  geom_point(aes(alpha = as.character(year)), shape = 3, size=5) +
  scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") + scale_colour_manual(name = "Cluster", values = cluster_palette) +
  scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
                                                                                                      "#5199ff", #john's island
                                                                                                      "#aac856", #skeleton lake
                                                                                                      "#a67e48" #snowgoose bay
  ), name = "Site") +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) + theme(axis.title = element_blank())
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_fun_cluster_rtsne2.svg"), plot=image,
       units="mm", width=300, height=300)


#all OTUs included in the functional mapping, DPCoA distance
func_otus <- read.table("D:/VirtualBox/VirtualBox Share/16S/Hazen16S/hazen_func_table_groups.txt", header = T, row.names = 1)
func_otu_names <- rowSums(func_otus)
func_otu_names <- names(which(func_otu_names > 0))
func_otu_only <- subset_taxa(phyloseqs[[i]], taxa_names(phyloseqs[[i]]) %in% func_otu_names)
func_otu_only <- prune_taxa(taxa_sums(func_otu_only) > 0, func_otu_only)
#midpoint root the tree before running DPCoA
phy_tree(func_otu_only) = midpoint(phy_tree(func_otu_only))
func_otu_only_diss <- DPCoA(func_otu_only)

set.seed(42)
rtsne_func_otu_only <- data.frame(Rtsne(func_otu_only_diss$RaoDis, is_distance = T, perplexity = 5)$Y)
rownames(rtsne_func_otu_only) <- sample_names(func_otu_only)
rtsne_func_otu_only_clusters <- hdbscan(rtsne_func_otu_only, xdist = vegdist(rtsne_func_otu_only, "euclidean"), minPts = 3)
rtsne_func_otu_only <- cbind(rtsne_func_otu_only, cluster = rtsne_func_otu_only_clusters$cluster, site = sample_data(func_otu_only)$Site, year = sample_data(func_otu_only)$Year, depth = sample_data(func_otu_only)$Depth.cm)
image <- ggplot(data = rtsne_func_otu_only, aes(x = X1, y = X2)) + stat_ellipse(level = 0.68, aes(color = factor(cluster)), linetype = 2, size = 1.5) +
  geom_point(aes(color = factor(cluster), fill = site, size = depth), shape = 21, stroke = 2) +
  geom_point(aes(alpha = as.character(year)), shape = 3, size=5) +
  scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") + scale_colour_manual(name = "Cluster", values = c("black", cluster_palette)) +
  scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
                                                                                                      "#5199ff", #john's island
                                                                                                      "#aac856", #skeleton lake
                                                                                                      "#a67e48" #snowgoose bay
  ), name = "Site") +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) + theme(axis.title = element_blank())
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_func_otu_cluster_rtsne2.svg"), plot=image,
       units="mm", width=300, height=300)

#check functional group differences between groups

#add lists to populate with random forest data
models <- list(NULL)
models_2 <- list(NULL)
pd_data <- list(NULL)
pd_data_2 <- list(NULL)
plot_pd_data <- list(NULL)
plot_pd_data_2 <- list(NULL)
tax_level_nums <- NA
full_model_errors <- list(NULL)
model_errors_2 <- list(NULL)
all_full_model_errors <- list(NULL)
all_model_errors_2 <- list(NULL)

#what are differences (both functional and taxonomic) between the functionally clustered groups
for (l in 1:2) {
  if (l == 1){predictors <- t(otu_table(func_pruned))} else {predictors <- t(otu_table(func_otu_only))}

  response <- factor(rtsne_func$cluster)

  rf_data <- data.frame(response, predictors)

  #run random forest model
  set.seed(42)
  classify <- ranger(response ~ ., data = rf_data, num.trees=5000, importance="impurity")


  #calculate model Kappa for the full model (inherently imbalanced data sets so we are using Cohen's Kappa to compare the models)
  pred1 <- classify$predictions
  kappa1 <- postResample(pred1, rf_data$response)[[2]]

  #sort the data by feature importance
  importances <- sort(importance(classify), decreasing = T)
  #restrict the number of variables to something reasonable to help with memory management
  if (length(importances) >= 250) {importances <- importances[1:250]}

  #add one feature at a time to the model in order of importance and compare kappas

  #reorder the sets
  rf_data_2 <- rf_data[c("response", names(importances))]

  comp_classify <- list(NULL)
  kappa2 <- NA
  for (k in 1:length(importances)) {
    #if all the importances are below the mean (all are 0?) break the loop
    if (k == 0) {break}
    new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
    set.seed(42)
    comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
    pred2 <- comp_classify[[k]]$predictions
    kappa2[k] <- postResample(pred2, new_data$response)[[2]]
    if (kappa2[k] == 1) {
      break
    }
  }
  if (length(kappa2) == length(importances)) {
    k <- min(which(kappa2 %in% max(kappa2)))
  }
  #this will store either values of the first "1" or the highest kappa value of all the models run
  all_model_errors_2[[j]] <- kappa2[k]


  #make partial dependence plots of the best model using all of the data
  nfeatures <- comp_classify[[k]]$num.independent.variables
  rf_data2 <- rf_data_2[1:(nfeatures+1)]
  pd <- partial_dependence(comp_classify[[k]], vars=colnames(rf_data2)[-1], data=rf_data2, n=c(25,nrow(rf_data2)))
  plot_pd_data <- plot_pd(pd)$data


  plot_pd_data$variable <- gsub("\\.", " ", plot_pd_data$variable)
  if (l == 1){taxonomy_glom <- prune_taxa(levels(plot_pd_data$variable), func_pruned)} else {taxonomy_glom <- prune_taxa(levels(plot_pd_data$variable), func_otu_only)}
  taxonomy_data.frame <- data.frame(variable=rownames(otu_table(taxonomy_glom)))
  taxonomy_merged <- merge(plot_pd_data,taxonomy_data.frame,by="variable")

  pd_ggplot <- ggplot(data = taxonomy_merged, aes(value, prediction*100)) + geom_line(aes(colour=variable), size= 1) +
    scale_x_continuous(trans="log2") + labs(x="Normalized abundance", y="Prediction (% chance to be classified)", colour="Function") +
    facet_grid(class~variable, scales="free") + theme(legend.position="none")

  if (l == 1) {
    ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/func_clustering_pd.svg"), plot=pd_ggplot,
           units="mm", width=300, height=200)
  } else  {
    ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/tax_clustering_pd.svg"), plot=pd_ggplot,
           units="mm", width=300, height=200)
  }
}

rm(comp_classify)
gc()

#check differences in physicochemistry between groups

#what are differences in physicochemistry of the samples between the clustered groups
for (l in 1:3) {
  if (l == 1){
    response <- factor(rtsne_common$cluster)
    cluster_data_name <- "tax"
  } else if (l == 2) {
    response <- factor(rtsne_func$cluster)
    cluster_data_name <- "func"
  } else if (l == 3) {
    response <- factor(rtsne_func_otu_only$cluster)
    cluster_data_name <- "func_otu"
  }

  predictors <- data.frame(sample_data(func_pruned))[-c(4,11)]
  predictors$core <- gsub("(.+)\\s[0-9].*","\\1",rownames(predictors))

  rf_data <- data.frame(response, predictors)

  #remove outlier data points
  if (any(levels(response) %in% 0)) {
    rf_data <- rf_data[-which(rf_data$response %in% 0),]
    rf_data$response <- factor(rf_data$response)
  }
  
  #run random forest model
  set.seed(42)
  classify <- ranger(response ~ ., data = rf_data, num.trees=5000, importance="impurity")


  #calculate model Kappa for the full model (inherently imbalanced data sets so we are using Cohen's Kappa to compare the models)
  pred1 <- classify$predictions
  kappa1 <- postResample(pred1, rf_data$response)[[2]]

  #sort the data by feature importance
  importances <- sort(importance(classify), decreasing = T)

  #add one feature at a time to the model in order of importance and compare kappas

  #reorder the sets
  rf_data_2 <- rf_data[c("response", names(importances))]

  comp_classify <- list(NULL)
  kappa2 <- NA
  for (k in 1:length(importances)) {
    #if all the importances are below the mean (all are 0?) break the loop
    if (k == 0) {break}
    new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
    set.seed(42)
    comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
    pred2 <- comp_classify[[k]]$predictions
    kappa2[k] <- postResample(pred2, new_data$response)[[2]]
    if (kappa2[k] == 1) {
      break
    }
  }
  if (length(kappa2) == length(importances)) {
    k <- min(which(kappa2 %in% max(kappa2)))
  }
  #this will store either values of the first "1" or the highest kappa value of all the models run
  all_model_errors_2[[j]] <- kappa2[k]


  #make partial dependence plots of the best model using all of the data
  nfeatures <- comp_classify[[k]]$num.independent.variables
  rf_data2 <- rf_data_2[1:(nfeatures+1)]
  pd <- partial_dependence(comp_classify[[k]], vars=colnames(rf_data2)[-1], data=rf_data2, n=c(25,nrow(rf_data2)))
  plot_pd_continuous <- plot_pd(pd)$data
  plot_pd_continuous <- na.omit(plot_pd_continuous)

  if (nrow(plot_pd_continuous) > 0) {
    #inset for the tax clustering
    if (l == 1) {
  plot_pd_continuous$variable <- gsub("\\.", " ", plot_pd_continuous$variable)
  rf_data2 <- cbind(rf_data2, sample_data(common)[,c(2,3)])
  rf_data2 <- rename(rf_data2, c("response" = "class"))
  plot_pd_continuous2 <- merge(plot_pd_continuous, rf_data2, by = "class")
  pd_continuous_ggplot2 <- ggplot(data = plot_pd_continuous2, aes(value, prediction*100)) + geom_line(aes(colour=class), size= 1.5) + #geom_vline(aes(xintercept=pH, colour=Site, linetype=Year), size = 1) +
    labs(x="value", y="Prediction (% chance to be classified)") + scale_color_manual(name = "Cluster", values = c(cluster_palette[1:length(levels(plot_pd_continuous$class))]
                                                                                                                 # , "#910a80", #deep hole
                                                                                                                 #  "#5199ff", #john's island
                                                                                                                 #  "#aac856", #skeleton lake
                                                                                                                 #  "#a67e48" #snowgoose bay
                                                                                                                  )) + 
    scale_linetype_manual(values=c(2,1)) + theme(legend.position="none")
    } 
    plot_pd_continuous$variable <- gsub("\\.", " ", plot_pd_continuous$variable)
    pd_continuous_ggplot <- ggplot(data = plot_pd_continuous, aes(value, prediction*100)) + geom_line(colour="black", size= 1) +
      labs(x="value", y="Prediction (% chance to be classified)") + facet_grid(class~variable, scales="free") + 
      theme(legend.position="none")
    plotheight  <- 200
  }

  #gather, plot and save categorical variable partial dependency plot (if applicable)
  if (any(colnames(pd) %in% c("core", "Site"))) {
  pd_categorical <- pd[, which(colnames(pd) %in% c("core", "Site", levels(response)))]
  pd_categorical <- pd_categorical[!apply(pd_categorical, 1, function (x) sum(is.na(x)) >= 2),]
  pd_categorical <- melt(pd_categorical)
  pd_categorical$type <- ifelse(is.na(pd_categorical$core),"site", "core")
  if (length(levels(factor(pd_categorical$type))) > 1) {
    pd_categorical <- as.data.frame(rbindlist(list(pd_categorical[which(is.na(pd_categorical$core)),-1],pd_categorical[which(is.na(pd_categorical$Site)),-2])))

    pd_categorical_ggplot <- ggplot(data = pd_categorical, aes(x = Site, y = value*100)) + geom_bar(stat="identity") +
      labs(y="Prediction (% chance to be classified", x="") + facet_grid(variable~type, scales = "free_x") +
      theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
  } else {
    pd_categorical_ggplot <- ggplot(data = pd_categorical, aes(x = core, y = value*100)) + geom_bar(stat="identity") +
      labs(y="Prediction (% chance to be classified", x="") + facet_grid(variable~., scales = "free_x") +
      theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
  }
  }


  if (nrow(plot_pd_continuous) > 0) {
    ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", cluster_data_name, "_clustering_pd_continuous.svg"), plot=pd_continuous_ggplot,
           units="mm", width=300, height=200)
    if (l == 1) {
      ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", cluster_data_name, "_clustering_pd_continuous_inset.svg"), plot=pd_continuous_ggplot2,
             units="mm", width=300, height=150)
    }
  }
  if (any(colnames(pd) %in% c("core", "Site"))) {
    ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", cluster_data_name, "_clustering_pd_categorical.svg"), plot=pd_categorical_ggplot,
           units="mm", width=300, height=200)
  }
}

}

# if (i == 1) {
# #clustering analysis on Lake Hazen samples only
# cluster_palette <- c("#00d092",
#                      "#d660cd",
#                      "#ae9300",
#                      "#0167bc",
#                      "#e65141",
#                      "#6c0034",
#                      "red",
#                      "green",
#                      "blue")
# 
# #all common >0.01% OTUs and DPCoA distance
# common_hazen_only <- subset_samples(common, Lake %in% "Lake Hazen")
# common_hazen_only <- prune_taxa(taxa_sums(common_hazen_only) > 0, common_hazen_only)
# #midpoint root the tree before running DPCoA
# phy_tree(common_hazen_only) = midpoint(phy_tree(common_hazen_only))
# common_hazen_only_diss <- DPCoA(common_hazen_only)
# common_hazen_only_matrix <- wisconsin(sqrt(veganifyOTU(common_hazen_only)))
# bestk <- silcheck(common_hazen_only_diss$RaoDis, diss=TRUE)[1]
# bestk2 <- hopach(data=common_hazen_only_matrix, dmat=common_hazen_only_diss$RaoDis, d= "euclid", clusters= "greedy")
# part <- pam(common_hazen_only_diss$RaoDis, k=bestk, pamonce = 2)
# part
# 
# sil_common_hazen_only <- data.frame(summary(part)$silinfo$widths)
# sil_common_hazen_only$sample <- rownames(sil_common_hazen_only)
# sil_common_hazen_only$sample <- factor(sil_common_hazen_only$sample, levels=sil_common_hazen_only$sample)
# sil_common_hazen_only_avg <- ddply(sil_common_hazen_only,~cluster,summarise,avg=paste0("avg. width = ", round(mean(sil_width), 2)))
# sil_common_hazen_only_avg <- merge(sil_common_hazen_only_avg, ddply(sil_common_hazen_only,~cluster,summarise,n=length(cluster)), by = "cluster")
# 
# image <- ggplot(data = sil_common_hazen_only, aes(x = sample, y = sil_width)) + geom_bar(stat = "identity") + facet_grid(cluster~., scales = "free_y") + coord_flip() + 
#   geom_text(data = sil_common_hazen_only_avg, aes(x = n, y = 0.66, label = avg), colour = "black", size = 3) + ylab("Silhouette width") + 
#   theme(axis.text.y = element_text(size=12))
# ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_common_silhouettes.svg"), plot=image,
#        units="mm", width=300, height=300)
# 
# pca_common_hazen_only <- prcomp(common_hazen_only_diss$RaoDis)
# image <- autoplot(pca_common_hazen_only, var.axes = F, alpha = 0, groups = factor(part$clustering)) +
#  stat_ellipse(level = 0.68, aes(color = factor(part$clustering)), linetype = 2, size = 1.5) +
#  geom_point(aes(color = factor(part$clustering), fill = sample_data(common_hazen_only)$Site, size = sample_data(common_hazen_only)$Depth.cm), shape = 21, stroke = 2) +
#  geom_point(aes(alpha = as.character(sample_data(common_hazen_only)$Year)), shape = 3, size=5) +
#  scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") + scale_colour_manual(name = "Cluster", values = cluster_palette) +
#  scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
#                                                                                                      "#5199ff", #john's island
#                                                                                                      "#a67e48" #snowgoose bay
#  ), name = "Site") +
#  guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) +
#  xlab(paste0("PC1 (", round(100*summary(pca_common_hazen_only)$importance[2],1), "% explained var.)")) + ylab(paste0("PC2 (", round(100*summary(pca_common_hazen_only)$importance[5],1), "% explained var.)"))
# 
# ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_common_cluster_pca.svg"), plot=image,
#       units="mm", width=300, height=300)
# 
# set.seed(42)
# rtsne_common_hazen_only <- data.frame(Rtsne(common_hazen_only_diss$RaoDis, is_distance = T, perplexity = 5)$Y)
# rownames(rtsne_common_hazen_only) <- sample_names(common_hazen_only)
# rtsne_common_hazen_only_clusters <- hdbscan(rtsne_common_hazen_only, xdist = vegdist(rtsne_common_hazen_only, "euclidean"), minPts = 3)
# rtsne_common_hazen_only <- cbind(rtsne_common_hazen_only, cluster = rtsne_common_hazen_only_clusters$cluster, site = sample_data(common_hazen_only)$Site, year = sample_data(common_hazen_only)$Year, depth = sample_data(common_hazen_only)$Depth.cm)
# image <- ggplot(data = rtsne_common_hazen_only, aes(x = X1, y = X2)) + stat_ellipse(level = 0.68, aes(color = factor(cluster)), linetype = 2, size = 1.5) + 
#    geom_point(aes(color = factor(cluster), fill = site, size = depth), shape = 21, stroke = 2) +
#    geom_point(aes(alpha = as.character(year)), shape = 3, size=5) +
#    scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") + scale_colour_manual(name = "Cluster", values = cluster_palette) +
#    scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
#                                                                                                        "#5199ff", #john's island
#                                                                                                        "#a67e48" #snowgoose bay
#    ), name = "Site") +
#    guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) + theme(axis.title = element_blank())
# ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_common_cluster_rtsne.svg"), plot=image,
#       units="mm", width=300, height=300)
# 
# 
# #all functional mapping groups, Bray-Curtis dissimilarity
# func_hazen_only <- subset_samples(func_pruned, Lake %in% "Lake Hazen")
# func_hazen_only <- prune_taxa(taxa_sums(func_hazen_only) > 0, func_hazen_only)
# func_hazen_only_matrix <- wisconsin(sqrt(veganifyOTU(func_hazen_only)))
# func_hazen_only_diss <- vegdist(func_hazen_only_matrix, distance = "bray")
# bestk_fun <- silcheck(func_hazen_only_diss, diss=TRUE)[1]
# bestk_fun2 <- hopach(data=func_hazen_only_matrix, dmat=func_hazen_only_diss, d= "euclid", clusters= "greedy")
# part_fun <- pam(func_hazen_only_diss, k=bestk_fun, pamonce = 2)
# part_fun
# 
# sil_fun_hazen_only <- data.frame(summary(part_fun)$silinfo$widths)
# sil_fun_hazen_only$sample <- rownames(sil_fun_hazen_only)
# sil_fun_hazen_only$sample <- factor(sil_fun_hazen_only$sample, levels=sil_fun_hazen_only$sample)
# sil_fun_hazen_only_avg <- ddply(sil_fun_hazen_only,~cluster,summarise,avg=paste0("avg. width = ", round(mean(sil_width), 2)))
# sil_fun_hazen_only_avg <- merge(sil_fun_hazen_only_avg, ddply(sil_fun_hazen_only,~cluster,summarise,n=length(cluster)), by = "cluster")
# 
# image <- ggplot(data = sil_fun_hazen_only, aes(x = sample, y = sil_width)) + geom_bar(stat = "identity") + facet_grid(cluster~., scales = "free_y") + coord_flip() + 
#   geom_text(data = sil_fun_hazen_only_avg, aes(x = n, y = 0.59, label = avg), colour = "black", size = 3) + ylab("Silhouette width") + 
#   theme(axis.text.y = element_text(size=12)) 
# ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_fun_silhouettes.svg"), plot=image,
#        units="mm", width=300, height=300)
# 
# pca_func_hazen_only <- prcomp(func_hazen_only_diss)
# image <- autoplot(pca_func_hazen_only, var.axes = F, alpha = 0, groups = factor(part_fun$clustering)) +
#   stat_ellipse(level = 0.68, aes(color = factor(part_fun$clustering)), linetype = 2, size = 1.5) +
#   geom_point(aes(color = factor(part_fun$clustering), fill = sample_data(func_hazen_only)$Site, size = sample_data(func_hazen_only)$Depth.cm), shape = 21, stroke = 2) +
#   geom_point(aes(alpha = as.character(sample_data(func_hazen_only)$Year)), shape = 3, size=4) +
#   scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") + scale_colour_manual(name = "Cluster", values = cluster_palette) +
#   scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
#                                                                                                       "#5199ff", #john's island
#                                                                                                       "#a67e48" #snowgoose bay
#   ), name = "Site") +
#   guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) +
#   xlab(paste0("PC1 (", round(100*summary(pca_func_hazen_only)$importance[2],1), "% explained var.)")) + ylab(paste0("PC2 (", round(100*summary(pca_func_hazen_only)$importance[5],1), "% explained var.)"))
# 
# ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_fun_cluster_pca.svg"), plot=image,
#        units="mm", width=300, height=300)
# 
# set.seed(42)
# rtsne_func_hazen_only <- data.frame(Rtsne(func_hazen_only_diss, is_distance = T, perplexity = 5)$Y)
# rtsne_func_hazen_only_clusters <- hdbscan(rtsne_func_hazen_only, xdist = vegdist(rtsne_func_hazen_only, "euclidean"), minPts = 3)
# rtsne_func_hazen_only <- cbind(rtsne_func_hazen_only, cluster = rtsne_func_hazen_only_clusters$cluster, site = sample_data(func_hazen_only)$Site, year = sample_data(func_hazen_only)$Year, depth = sample_data(func_hazen_only)$Depth.cm)
# image <- ggplot(data = rtsne_func_hazen_only, aes(x = X1, y = X2)) + stat_ellipse(level = 0.68, aes(color = factor(cluster)), linetype = 2, size = 1.5) + 
#   geom_point(aes(color = factor(cluster), fill = site, size = depth), shape = 21, stroke = 2) +
#   geom_point(aes(alpha = as.character(year)), shape = 3, size=5) +
#   scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") + scale_colour_manual(name = "Cluster", values = c("black", cluster_palette)) +
#   scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
#                                                                                                       "#5199ff", #john's island
#                                                                                                       "#a67e48" #snowgoose bay
#   ), name = "Site") +
#   guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) + theme(axis.title = element_blank())
# ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_fun_cluster_rtsne.svg"), plot=image,
#        units="mm", width=300, height=300)
# 
# 
# #all OTUs included in the functional mapping, DPCoA distance
# func_otus <- read.table("D:/VirtualBox/VirtualBox Share/16S/Hazen16S/hazen_func_table_groups.txt", header = T, row.names = 1)
# func_otu_names <- rowSums(func_otus)
# func_otu_names <- names(which(func_otu_names > 0))
# func_otu_only <- subset_taxa(phyloseqs[[i]], taxa_names(phyloseqs[[i]]) %in% func_otu_names)
# func_otu_only <- subset_samples(func_otu_only, Lake %in% "Lake Hazen")
# 
# func_otu_only <- prune_taxa(taxa_sums(func_otu_only) > 0, func_otu_only)
# #midpoint root the tree before running DPCoA
# phy_tree(func_otu_only) = midpoint(phy_tree(func_otu_only))
# func_otu_only_diss <- DPCoA(func_otu_only)
# func_otu_only_matrix <- wisconsin(sqrt(veganifyOTU(func_otu_only)))
# bestk_func_otu <- silcheck(func_otu_only_diss$RaoDis, diss=TRUE)[1]
# bestk2_func_otu <- hopach(data=func_otu_only_matrix, dmat=func_otu_only_diss$RaoDis, d= "euclid", clusters= "greedy")
# part_func_otu <- pam(func_otu_only_diss$RaoDis, k=bestk_func_otu, pamonce = 2)
# part_func_otu
# 
# sil_func_otu_only <- data.frame(summary(part_func_otu)$silinfo$widths)
# sil_func_otu_only$sample <- rownames(sil_func_otu_only)
# sil_func_otu_only$sample <- factor(sil_func_otu_only$sample, levels=sil_func_otu_only$sample)
# sil_func_otu_only_avg <- ddply(sil_func_otu_only,~cluster,summarise,avg=paste0("avg. width = ", round(mean(sil_width), 2)))
# sil_func_otu_only_avg <- merge(sil_func_otu_only_avg, ddply(sil_func_otu_only,~cluster,summarise,n=length(cluster)), by = "cluster")
# 
# image <- ggplot(data = sil_func_otu_only, aes(x = sample, y = sil_width)) + geom_bar(stat = "identity") + facet_grid(cluster~., scales = "free_y") + coord_flip() + 
#   geom_text(data = sil_func_otu_only_avg, aes(x = n, y = 0.6, label = avg), colour = "black", size = 3) + ylab("Silhouette width") + 
#   theme(axis.text.y = element_text(size=12))
# ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_func_otu_silhouettes.svg"), plot=image,
#        units="mm", width=300, height=300)
# 
# pca_func_otu_only <- prcomp(func_otu_only_diss$RaoDis)
# image <- autoplot(pca_func_otu_only, var.axes = F, alpha = 0, groups = factor(part_func_otu$clustering)) +
#   stat_ellipse(level = 0.68, aes(color = factor(part$clustering)), linetype = 2, size = 1.5) +
#   geom_point(aes(color = factor(part$clustering), fill = sample_data(func_otu_only)$Site, size = sample_data(func_otu_only)$Depth.cm), shape = 21, stroke = 2) +
#   geom_point(aes(alpha = as.character(sample_data(func_otu_only)$Year)), shape = 3, size=4) +
#   scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") + scale_colour_manual(name = "Cluster", values = cluster_palette) +
#   scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
#                                                                                                       "#5199ff", #john's island
#                                                                                                       "#a67e48" #snowgoose bay
#   ), name = "Site") +
#   guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) +
#   xlab(paste0("PC1 (", round(100*summary(pca_common_hazen_only)$importance[2],1), "% explained var.)")) + ylab(paste0("PC2 (", round(100*summary(pca_common_hazen_only)$importance[5],1), "% explained var.)"))
# 
# ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_func_otu_cluster_pca.svg"), plot=image,
#        units="mm", width=300, height=300)
# 
# set.seed(42)
# rtsne_func_otu_only <- data.frame(Rtsne(func_otu_only_diss$RaoDis, is_distance = T, perplexity = 5)$Y)
# rownames(rtsne_func_otu_only) <- sample_names(func_otu_only)
# rtsne_func_otu_only_clusters <- hdbscan(rtsne_func_otu_only, xdist = vegdist(rtsne_func_otu_only, "euclidean"), minPts = 3)
# rtsne_func_otu_only <- cbind(rtsne_func_otu_only, cluster = rtsne_func_otu_only_clusters$cluster, site = sample_data(func_otu_only)$Site, year = sample_data(func_otu_only)$Year, depth = sample_data(func_otu_only)$Depth.cm)
# image <- ggplot(data = rtsne_func_otu_only, aes(x = X1, y = X2)) + stat_ellipse(level = 0.68, aes(color = factor(cluster)), linetype = 2, size = 1.5) + 
#   geom_point(aes(color = factor(cluster), fill = site, size = depth), shape = 21, stroke = 2) +
#   geom_point(aes(alpha = as.character(year)), shape = 3, size=5) +
#   scale_alpha_manual(values = c("2014" = 1, "2015" = 0), name = "Year") + scale_colour_manual(name = "Cluster", values = cluster_palette) +
#   scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
#                                                                                                       "#5199ff", #john's island
#                                                                                                       "#a67e48" #snowgoose bay
#   ), name = "Site") +
#   guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) + theme(axis.title = element_blank())
# ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_func_otu_cluster_rtsne.svg"), plot=image,
#        units="mm", width=300, height=300)
# 
# #check functional group differences between groups
# 
# #add lists to populate with random forest data
# models <- list(NULL)
# models_2 <- list(NULL)
# pd_data <- list(NULL)
# pd_data_2 <- list(NULL)
# plot_pd_data <- list(NULL)
# plot_pd_data_2 <- list(NULL)
# tax_level_nums <- NA
# full_model_errors <- list(NULL)
# model_errors_2 <- list(NULL)
# all_full_model_errors <- list(NULL)
# all_model_errors_2 <- list(NULL)
# 
# #what are differences (both functional and taxonomic) between the functionally clustered groups
# for (l in 1:2) {
#   if (l == 1){predictors <- t(otu_table(func_hazen_only))} else {predictors <- t(otu_table(func_otu_only))}
#   
#   response <- factor(rtsne_func_hazen_only$cluster)
#   
#   rf_data <- data.frame(response, predictors)
#   
#   #run random forest model 
#   set.seed(42)
#   classify <- ranger(response ~ ., data = rf_data, num.trees=5000, importance="impurity")
#   
#   
#   #calculate model Kappa for the full model (inherently imbalanced data sets so we are using Cohen's Kappa to compare the models)
#   pred1 <- classify$predictions
#   kappa1 <- postResample(pred1, rf_data$response)[[2]]
#   
#   #sort the data by feature importance
#   importances <- sort(importance(classify), decreasing = T)
#   #restrict the number of variables to something reasonable to help with memory management
#   if (length(importances) >= 250) {importances <- importances[1:250]}
#   
#   #add one feature at a time to the model in order of importance and compare kappas
#   
#   #reorder the sets
#   rf_data_2 <- rf_data[c("response", names(importances))]
#   
#   comp_classify <- list(NULL)
#   kappa2 <- NA
#   for (k in 1:length(importances)) {
#     #if all the importances are below the mean (all are 0?) break the loop
#     if (k == 0) {break}
#     new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
#     set.seed(42)
#     comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
#     pred2 <- comp_classify[[k]]$predictions
#     kappa2[k] <- postResample(pred2, new_data$response)[[2]]
#     if (kappa2[k] == 1) {
#       break
#     }
#   }
#   if (length(kappa2) == length(importances)) {
#     k <- min(which(kappa2 %in% max(kappa2)))
#   }
#   #this will store either values of the first "1" or the highest kappa value of all the models run
#   all_model_errors_2[[j]] <- kappa2[k]
#   
#   
#   #make partial dependence plots of the best model using all of the data
#   nfeatures <- comp_classify[[k]]$num.independent.variables
#   rf_data2 <- rf_data_2[1:(nfeatures+1)]
#   pd <- partial_dependence(comp_classify[[k]], vars=colnames(rf_data2)[-1], data=rf_data2, n=c(25,nrow(rf_data2)))
#   plot_pd_data <- plot_pd(pd)$data
#   
#   
#   plot_pd_data$variable <- gsub("\\.", " ", plot_pd_data$variable)
#   if (l == 1){taxonomy_glom <- prune_taxa(levels(plot_pd_data$variable), func_hazen_only)} else {taxonomy_glom <- prune_taxa(levels(plot_pd_data$variable), func_otu_only)}
#   taxonomy_data.frame <- data.frame(variable=rownames(otu_table(taxonomy_glom)))
#   taxonomy_merged <- merge(plot_pd_data,taxonomy_data.frame,by="variable")
#   
#   pd_ggplot <- ggplot(data = taxonomy_merged, aes(value, prediction*100)) + geom_line(aes(colour=variable), size= 1) +
#     scale_x_continuous(trans="log2") + labs(x="Normalized abundance", y="Prediction (% chance to be classified)", colour="Function") +
#     facet_grid(class~variable, scales="free") + theme(legend.position="none")
#   
#   if (l == 1) {
#     ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/func_clustering_pd.svg"), plot=pd_ggplot,
#            units="mm", width=1200, height=200)
#   } else  {
#     ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/tax_clustering_pd.svg"), plot=pd_ggplot,
#            units="mm", width=1200, height=200)
#   }
# }
# 
# rm(comp_classify)
# gc()
# 
# #check differences in physicochemistry between groups
# 
# #what are differences in physicochemistry of the samples between the clustered groups
# for (l in 1:3) {
#   if (l == 1){
#     response <- factor(rtsne_common_hazen_only$cluster)
#     cluster_data_name <- "tax"
#   } else if (l == 2) {
#     response <- factor(rtsne_func_hazen_only$cluster)
#     cluster_data_name <- "func"
#   } else if (l == 3) {
#     response <- factor(rtsne_func_otu_only$cluster)
#     cluster_data_name <- "func_otu"
#   }
#   
#   predictors <- data.frame(sample_data(func_hazen_only))[-c(5,6,7)]
#   predictors$core <- gsub("(.+)\\s[0-9].*","\\1",rownames(predictors))
#   
#   rf_data <- data.frame(response, predictors)
#   
#   #run random forest model 
#   set.seed(42)
#   classify <- ranger(response ~ ., data = rf_data, num.trees=5000, importance="impurity")
#   
#   
#   #calculate model Kappa for the full model (inherently imbalanced data sets so we are using Cohen's Kappa to compare the models)
#   pred1 <- classify$predictions
#   kappa1 <- postResample(pred1, rf_data$response)[[2]]
#   
#   #sort the data by feature importance
#   importances <- sort(importance(classify), decreasing = T)
#   
#   #add one feature at a time to the model in order of importance and compare kappas
#   
#   #reorder the sets
#   rf_data_2 <- rf_data[c("response", names(importances))]
#   
#   comp_classify <- list(NULL)
#   kappa2 <- NA
#   for (k in 1:length(importances)) {
#     #if all the importances are below the mean (all are 0?) break the loop
#     if (k == 0) {break}
#     new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
#     set.seed(42)
#     comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
#     pred2 <- comp_classify[[k]]$predictions
#     kappa2[k] <- postResample(pred2, new_data$response)[[2]]
#     if (kappa2[k] == 1) {
#       break
#     }
#   }
#   if (length(kappa2) == length(importances)) {
#     k <- min(which(kappa2 %in% max(kappa2)))
#   }
#   #this will store either values of the first "1" or the highest kappa value of all the models run
#   all_model_errors_2[[j]] <- kappa2[k]
#   
#   
#   #make partial dependence plots of the best model using all of the data
#   nfeatures <- comp_classify[[k]]$num.independent.variables
#   rf_data2 <- rf_data_2[1:(nfeatures+1)]
#   pd <- partial_dependence(comp_classify[[k]], vars=colnames(rf_data2)[-1], data=rf_data2, n=c(25,nrow(rf_data2)))
#   plot_pd_continuous <- plot_pd(pd)$data
#   plot_pd_continuous <- na.omit(plot_pd_continuous)
#   
#   if (nrow(plot_pd_continuous) > 0) {
#   plot_pd_continuous$variable <- gsub("\\.", " ", plot_pd_continuous$variable)
#   pd_continuous_ggplot <- ggplot(data = plot_pd_continuous, aes(value, prediction*100)) + geom_line(aes(colour="black"), size= 1) +
#     labs(x="value", y="Prediction (% chance to be classified)", colour="Function") +
#     facet_grid(class~variable, scales="free") + theme(legend.position="none")
#   }
#   
#   #gather, plot and save categorical variable partial dependency plot (if applicable)
#   if (any(colnames(pd) %in% c("core", "Site"))) {
#   pd_categorical <- pd[, which(colnames(pd) %in% c("core", "Site", levels(response)))]
#   pd_categorical <- pd_categorical[!apply(pd_categorical, 1, function (x) sum(is.na(x)) >= 2),]
#   pd_categorical <- melt(pd_categorical)
#   pd_categorical$type <- ifelse(is.na(pd_categorical$core),"site", "core")
#   if (length(levels(factor(pd_categorical$type))) > 1) {
#     pd_categorical <- as.data.frame(rbindlist(list(pd_categorical[which(is.na(pd_categorical$core)),-1],pd_categorical[which(is.na(pd_categorical$Site)),-2])))
#   
#     pd_categorical_ggplot <- ggplot(data = pd_categorical, aes(x = Site, y = value*100)) + geom_bar(stat="identity") +
#       labs(y="Prediction (% chance to be classified", x="") + facet_grid(variable~type, scales = "free_x") + 
#       theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
#   } else {
#     pd_categorical_ggplot <- ggplot(data = pd_categorical, aes(x = core, y = value*100)) + geom_bar(stat="identity") +
#       labs(y="Prediction (% chance to be classified", x="") + facet_grid(variable~., scales = "free_x") + 
#       theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
#   }
#   }
#   
#   
#   if (nrow(plot_pd_continuous) > 0) {
#     ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", cluster_data_name, "_clustering_pd_continuous.svg"), plot=pd_continuous_ggplot,
#            units="mm", width=300, height=200)
#   } 
#   if (any(colnames(pd) %in% c("core", "Site"))) {
#     ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", cluster_data_name, "_clustering_pd_categorical.svg"), plot=pd_categorical_ggplot,
#            units="mm", width=300, height=200)
#   }
# }
# 
# }

#make tables of data available for models
datasets <- c(rep("common",ncol(sample_data(common))),rep("func_pruned",ncol(sample_data(common))))
variables <- rep(colnames(sample_data(common)), 2)
datasets_ord <- c(rep("ordu_common",ncol(sample_data(common))),rep("ordu_func",ncol(sample_data(common))))
modelstorun_frame <- data.frame(set=datasets, model=variables, dataset_ord=datasets_ord)

#readd "n.seq" to continuous variables
continuous_variables <- c(continuous_variables,"n.seq")

#run ordisurf for all the continuous variables and extract all pvalues for both data sets
ordisurfs <- list(NULL)
pvalues <- c(NULL)
for (j in 1:nrow(modelstorun_frame)) {
  if (modelstorun_frame$model[j] %in% continuous_variables) {
    ordisurfs[[j]] <- ordisurf(eval(parse(text = paste0(modelstorun_frame$dataset_ord[j],"~",modelstorun_frame$model[j]))), data = data.frame(eval(parse(text = paste0("sample_data(",modelstorun_frame$set[j],")")))),
                               method = "REML", select = TRUE, penalty= 1.4, plot = F, scaling = 3, w =NULL, permu = 10000)
    pvalues[j] <- summary(ordisurfs[[j]])$p.table[1,4]
  } else {
    ordisurfs[[j]] <- NULL
    if (modelstorun_frame$model[j] %in% "common") {
      pvalues[j] <- unname(ef_common$factors$pvals[names(ef_common$factors$pvals) %in% modelstorun_frame$model[j]])
    } else {
      pvalues[j] <- unname(ef_func$factors$pvals[names(ef_func$factors$pvals) %in% modelstorun_frame$model[j]])
    }
  }
}

#add physicochemical vector fit pvalues to the end of the pvalues
pvalues_physchem <- c(ef_common$vectors$pvals,ef_func$vectors$pvals)

#get the mantel test pvalues
pvalues_mantel <- c(common_mantel$signif, func_mantel$signif)

#correct p-values for all variables, separate them for ordisurfs/categories, vectors, and mantel tests, and add them to the modelstorun frame
n_bonferroni <- length(c(pvalues,pvalues_physchem,pvalues_mantel))  #set the bonferroni n tests to number of all significance tests run per data set (surface fits + vectors + group permutations + mantel tests)
pvalues_bonferroni_corrected <- p.adjust(c(pvalues,pvalues_physchem,pvalues_mantel), method = "bonferroni", n = n_bonferroni)
pvalues_physchem_corrected <- pvalues_bonferroni_corrected[(length(pvalues)+1):(length(pvalues_bonferroni_corrected)-2)]
pvalues_mantel_corrected <- pvalues_bonferroni_corrected[(length(pvalues_bonferroni_corrected)-1):length(pvalues_bonferroni_corrected)]
pvalues_bonferroni_corrected <- pvalues_bonferroni_corrected[1:length(pvalues)]
modelstorun_frame$bonf_pvalues <- pvalues_bonferroni_corrected
modelstorun_frame$bonf_vector_pvalues <- NA
modelstorun_frame[modelstorun_frame$model %in% continuous_variables[-which(continuous_variables %in% "n.seq")],]$bonf_vector_pvalues <- pvalues_physchem_corrected

#gather mantel test data
mantels[[i]] <- data.frame(dataset=dataset_names[i],
                           set=c("common", "func"),
                           r2=c(common_mantel$statistic,func_mantel$statistic),
                           bonf_pvalues=pvalues_mantel_corrected)


#make a frame for the vectors
physchem_vectors <- data.frame(set = c(rep("OTUs", length(pvalues_physchem)/2), rep("FAPROTAX", length(pvalues_physchem)/2)),
                               model = c(names(ef_common$vectors$r), names(ef_func$vectors$r)),
                               bonf_pvalues = pvalues_physchem_corrected)
physchem_vectors <- cbind(physchem_vectors,(rbind(scores(ef_common, display="vectors"), scores(ef_func, display="vectors"))))

#retain only vectors with < 0.05 p-values
physchem_vectors <- physchem_vectors[physchem_vectors$bonf_pvalues < 0.05,]

#agglomerate taxa to data sets until Order level and remove groups that have less than 0.1% overall abundance
tax_levels <- list(NULL)
for (h in 2:4) {
  tax_levels[[h]] <- tax_glom(phyloseqs[[i]],taxrank = rank_names(phyloseqs[[i]])[h])
  #remove the tree slot because it causes errors later on (construct a new object without the tree)
  tax_levels[[h]] <- phyloseq(otu_table(tax_levels[[h]]), sample_data(tax_levels[[h]]), tax_table(tax_levels[[h]]))
  }

plot_pd_data <- list(NULL)

for (j in 1:nrow(modelstorun_frame)) {
  if (modelstorun_frame$set[j] %in% "common"){ #this run for taxonomic data to choose the optimal model among all of the levels
    tax_models <- list(NULL)
    model_errors <- list(NULL)
    tax_models_2 <- list(NULL)
    model_errors_2 <- list(NULL)
    rf_datas <- list(NULL)
    for (q in 2:length(tax_levels)) {
      #prepare data for input to random forest
      predictors <-
        t(otu_table(tax_levels[[q]]))
      if (modelstorun_frame$model[j] %in% continuous_variables) {
        response <-
          eval(parse(
            text = paste0(
              "sample_data(",
              modelstorun_frame$set[j],
              ")$",
              modelstorun_frame$model[j]
            )
          ))
        
      } else {
        response <-
          as.factor(eval(parse(
            text = paste0(
              "sample_data(",
              modelstorun_frame$set[j],
              ")$",
              modelstorun_frame$model[j]
            )
          )))
      }
      rf_data <- data.frame(response, predictors)
      
      #run random forest model 
      set.seed(42)
      classify <- ranger(response ~ ., data = rf_data, num.trees=5000, importance="impurity")
      tax_models[[q]] <- classify
      
      if (modelstorun_frame$model[j] %in% continuous_variables) {
        
        #calculate MSPE estimate with CI for the full model
        pred1 <- classify$predictions
        y <-  rf_data$response
        n <- length(y)
        
        # psi is the mean squared prediction error (MSPE) estimate
        # sigma2 is the estimate of the variance of the MSPE
        psi1 <- mean((y - pred1)^2)
        sigma21 <- 1/n * var((y - pred1)^2) 
        # 95% CI:
        full_MSPE <- c(psi1 - 1.96 * sqrt(sigma21), psi1, psi1 + 1.96 * sqrt(sigma21))
        
        #save errors for the full model
        full_model_errors[[q]] <- full_MSPE
        
        #sort the data by feature importance
        importances <- sort(importance(classify), decreasing = T)
        
        #add one feature at a time to the model in order of importance and calculate MSPEs for all of the models
        
        #reorder the sets
        rf_data_2 <- rf_data[c("response", names(importances))]
        rf_datas[[q]] <- rf_data_2
        
        comp_classify <- list(NULL)
        comp_MSPEs <- list(NULL)
        for (k in 1:length(importances)) {
          #if all the importances are below the mean (all are 0?) break the loop
          if (k == 0) {break}
          new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
          set.seed(42)
          comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
          pred2 <- comp_classify[[k]]$predictions
          y2 <- new_data$response
          n2 <- length(y2)
          psi2 <- mean((y2 - pred2)^2)
          sigma22 <- 1/n * var((y2 - pred2)^2) 
          # 95% CI:
          comp_MSPEs[[k]] <- c(psi2 - 1.96 * sqrt(sigma22), psi2, psi2 + 1.96 * sqrt(sigma22))
        }
        #find the minimum MSPE and get the feature names
        comp_MSPEs <- do.call("rbind", comp_MSPEs)
        nfeatures <- min(which(comp_MSPEs[,2] %in% min(comp_MSPEs[,2])))
        feature_names <- names(importances)[1:nfeatures]
        tax_models_2[[q]] <- comp_classify[[nfeatures]]
        model_errors_2[[q]] <- comp_MSPEs[nfeatures,]
      } else {
        #calculate model Kappa for the full model (inherently imbalanced data sets so we are using kappa to compare the models)
        pred1 <- classify$predictions
        kappa1 <- postResample(pred1, rf_data$response)[[2]]
        full_model_errors[[q]] <- kappa1
        
        #sort the data by feature importance
        importances <- sort(importance(classify), decreasing = T)
        
        #add one feature at a time to the model in order of importance and calculate Kappa values for all of the models
        
        #reorder the sets
        rf_data_2 <- rf_data[c("response", names(importances))]
        rf_datas[[q]] <- rf_data_2
        
        comp_classify <- list(NULL)
        kappa2 <- NA
        for (k in 1:length(importances)) {
          #if all the importances are below the mean (all are 0?) break the loop
          if (k == 0) {break}
          new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
          set.seed(42)
          comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
          pred2 <- comp_classify[[k]]$predictions
          kappa2[k] <- postResample(pred2, new_data$response)[[2]]
          if (kappa2[k] == 1) {
            break
          }
        }
        if (length(kappa2) == length(importances)) {
          k <- min(which(kappa2 %in% max(kappa2)))
        }
        #this will store either values of the first "1" or the highest kappa value of all the models run
        model_errors_2[[q]] <- kappa2[k]
        tax_models_2[[q]] <- comp_classify[[k]]
      }
    }
    
    if (modelstorun_frame$model[j] %in% continuous_variables) {
      tax_MSPEs <- unlist(sapply(model_errors_2, "[[", 2))
      min_tax_MSPE <- min(tax_MSPEs)
      tax_level_nums[j] <- min(which(tax_MSPEs %in% min_tax_MSPE)+1)
    } else {
      tax_kappas <- unlist(model_errors_2)
      max_tax_kappa <- max(tax_kappas)
      tax_level_nums[j] <- min(which(tax_kappas %in% max_tax_kappa)+1)
    }
    
    #save the best fit original model, new model, and lists of *all* model errors
    models[[j]] <- tax_models[[tax_level_nums[j]]]
    models_2[[j]] <-  tax_models_2[[tax_level_nums[j]]]
    all_full_model_errors[[j]] <- full_model_errors[[tax_level_nums[j]]]
    all_model_errors_2[[j]] <- model_errors_2[[tax_level_nums[j]]]

    #make partial dependence plots of the best model using all of the data
    nfeatures <- models_2[[j]]$num.independent.variables
    rf_data2 <- rf_datas[[tax_level_nums[j]]][1:(nfeatures+1)]
    pd <- partial_dependence(models_2[[j]], vars=colnames(rf_data2)[-1], data=rf_data2, n=c(25,nrow(rf_data2)))
    pd_data[[j]] <- pd
    plot_pd_data[[j]] <- plot_pd(pd)$data
    
  } else { #this run for functional data
    #prepare data for input to random forest
    predictors <-
      t(otu_table(eval(parse(
        text = paste0(modelstorun_frame$set[j])
      ))))
    if (modelstorun_frame$model[j] %in% continuous_variables) {
      response <-
        eval(parse(
          text = paste0(
            "sample_data(",
            modelstorun_frame$set[j],
            ")$",
            modelstorun_frame$model[j]
          )
        ))
      
    } else {
      response <-
        as.factor(eval(parse(
          text = paste0(
            "sample_data(",
            modelstorun_frame$set[j],
            ")$",
            modelstorun_frame$model[j]
          )
        )))
    }
    rf_data <- data.frame(response, predictors)
    
    #run random forest model 
    set.seed(42)
    classify <- ranger(response ~ ., data = rf_data, num.trees=5000, importance="impurity")
    models[[j]] <- classify
    
    if (modelstorun_frame$model[j] %in% continuous_variables) {
      
      #calculate MSPE estimate and CI for the full model
      pred1 <- classify$predictions
      y <-  rf_data$response
      n <- length(y)
      
      # psi is the mean squared prediction error (MSPE) estimate
      # sigma2 is the estimate of the variance of the MSPE
      psi1 <- mean((y - pred1)^2)
      sigma21 <- 1/n * var((y - pred1)^2) 
      # 95% CI:
      full_MSPE <- c(psi1 - 1.96 * sqrt(sigma21), psi1, psi1 + 1.96 * sqrt(sigma21))
      
      #save errors for the full model
      all_full_model_errors[[j]] <- full_MSPE
      
      #sort the data by feature importance
      importances <- sort(importance(classify), decreasing = T)
      
      #add one feature at a time to the model in order of importance and compare MSPEs (consider only features with > mean importance)
      #when a minimum in MSPE is reached, stop and return the model that was run 
      
      #reorder the sets
      rf_data_2 <- rf_data[c("response", names(importances))]
      
      comp_classify <- list(NULL)
      comp_MSPEs <- list(NULL)
      for (k in 1:length(importances)) {
        #if all the importances are below the mean (all are 0?) break the loop
        if (k == 0) {break}
        new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
        set.seed(42)
        comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
        pred2 <- comp_classify[[k]]$predictions
        y2 <- new_data$response
        n2 <- length(y2)
        psi2 <- mean((y2 - pred2)^2)
        sigma22 <- 1/n * var((y2 - pred2)^2) 
        # 95% CI:
        comp_MSPEs[[k]] <- c(psi2 - 1.96 * sqrt(sigma22), psi2, psi2 + 1.96 * sqrt(sigma22))
      }
      #find the minimum MSPE and get the feature names
      comp_MSPEs <- do.call("rbind", comp_MSPEs)
      nfeatures <- min(which(comp_MSPEs[,2] %in% min(comp_MSPEs[,2])))
      feature_names <- names(importances)[1:nfeatures]
      models_2[[j]] <- comp_classify[[nfeatures]]
      all_model_errors_2[[j]] <- comp_MSPEs[nfeatures,]
    } else {
      #calculate model Kappa for the full model (inherently imbalanced data sets so we are using Cohen's Kappa to compare the models)
      pred1 <- classify$predictions
      kappa1 <- postResample(pred1, rf_data$response)[[2]]
      all_full_model_errors[[j]] <- kappa1
      
      #sort the data by feature importance
      importances <- sort(importance(classify), decreasing = T)
      
      #add one feature at a time to the model in order of importance and compare Kappa values
      
      #reorder the sets
      rf_data_2 <- rf_data[c("response", names(importances))]
      
      comp_classify <- list(NULL)
      kappa2 <- NA
      for (k in 1:length(importances)) {
        #if all the importances are below the mean (all are 0?) break the loop
        if (k == 0) {break}
        new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
        set.seed(42)
        comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
        pred2 <- comp_classify[[k]]$predictions
        kappa2[k] <- postResample(pred2, new_data$response)[[2]]
        if (kappa2[k] == 1) {
          break
        }
      }
      if (length(kappa2) == length(importances)) {
        k <- min(which(kappa2 %in% max(kappa2)))
      }
      #this will store either values of the first "1" or the highest Kappa value of all the models run
      all_model_errors_2[[j]] <- kappa2[k]
      models_2[[j]] <- comp_classify[[k]]
      
    }
    
    #make partial dependence plots of the best model using all of the data
    nfeatures <- models_2[[j]]$num.independent.variables
    rf_data2 <- rf_data_2[1:(nfeatures+1)]
    pd <- partial_dependence(models_2[[j]], vars=colnames(rf_data2)[-1], data=rf_data2, n=c(25,nrow(rf_data2)))
    pd_data[[j]] <- pd
    plot_pd_data[[j]] <- plot_pd(pd)$data
  }
}

#format and save partial dependency ggplots
tax_level_sets <- colnames(tax)[1:4]
pd_ggplots <- list(NULL)
for (k in 1:nrow(modelstorun_frame)) {
  if (modelstorun_frame$set[k] %in% "common") {
    taxonomy_glom <- prune_taxa(levels(plot_pd_data[[k]]$variable), tax_levels[[tax_level_nums[k]]])
    #get a clearer taxonomy for the class levels
      taxonomy_cat <- list(NULL)
        for (l in 1:ntaxa(taxonomy_glom)){
        taxonomy_cat[[l]] <- tax_table(taxonomy_glom)[l,2:4]
        taxonomy_cat[[l]] <- gsub("NA","",taxonomy_cat[[l]])
        }
      taxonomy_data.frame <- data.frame(do.call("rbind", taxonomy_cat))
      taxonomy_data.frame$variable <- rownames(taxonomy_data.frame)

    taxonomy_merged <- merge(plot_pd_data[[k]],taxonomy_data.frame,by="variable")
    taxonomy_merged$taxonomy <- taxonomy_merged$taxonomy
    facet_groups <- capture.output(cat(tax_level_sets[2:tax_level_nums[k]],sep="+"))
    
    if (modelstorun_frame$model[k] %in% continuous_variables) {
    pd_ggplots[[k]] <- ggplot(data = taxonomy_merged, aes(value, response)) + geom_line(aes(colour=variable), size= 1) +
    scale_x_continuous(trans="log2") + labs(x="Normalized abundance", y=modelstorun_frame$model[k], colour=tax_level_sets[tax_level_nums[k]]) +
      theme(legend.position="none") 
    if (length(unique(taxonomy_merged$variable)) >= 30) {
      pd_ggplots[[k]] <- pd_ggplots[[k]] + eval(parse(text=paste0("facet_wrap(~",facet_groups,", scales=\"free_x\", labeller = label_context, ncol=16)")))
    } else if (length(unique(taxonomy_merged$variable)) >= 15) {
      pd_ggplots[[k]] <- pd_ggplots[[k]] + eval(parse(text=paste0("facet_wrap(~",facet_groups,", scales=\"free_x\", labeller = label_context, ncol=8)")))
    } else if (length(unique(taxonomy_merged$variable)) >= 5) {
      pd_ggplots[[k]] <- pd_ggplots[[k]] + eval(parse(text=paste0("facet_wrap(~",facet_groups,", scales=\"free_x\", labeller = label_context, ncol=4)")))
    } else {
      pd_ggplots[[k]] <- pd_ggplots[[k]] + eval(parse(text=paste0("facet_wrap(~",facet_groups,", scales=\"free_x\", labeller = label_context, ncol=2)")))
    }
    } else {
      pd_ggplots[[k]] <- ggplot(data = taxonomy_merged, aes(value, prediction*100)) + geom_line(aes(colour=variable), size= 1) +
        scale_x_continuous(trans="log2") + labs(x="Normalized abundance", y="Prediction (% chance to be classified)", colour=tax_level_sets[tax_level_nums[k]]) +
        theme(legend.position="none") + eval(parse(text=paste0("facet_grid(class~",facet_groups,", scales=\"free\", labeller = label_context)")))
    }
  } else {
    plot_pd_data[[k]]$variable <- gsub("\\.", " ", plot_pd_data[[k]]$variable)
    taxonomy_glom <- prune_taxa(levels(plot_pd_data[[k]]$variable), func_pruned)
    taxonomy_data.frame <- data.frame(variable=rownames(otu_table(taxonomy_glom)))
    taxonomy_merged <- merge(plot_pd_data[[k]],taxonomy_data.frame,by="variable")
    if (modelstorun_frame$model[k] %in% continuous_variables) {
      pd_ggplots[[k]] <- ggplot(data = taxonomy_merged, aes(value, response)) + geom_line(aes(colour=variable), size= 1) +
      scale_x_continuous(trans="log2") + labs(x="Normalized abundance", y=modelstorun_frame$model[k], colour="Function") +
        theme(legend.position="none")
      if (length(unique(taxonomy_merged$variable)) >= 15) {
        pd_ggplots[[k]] <- pd_ggplots[[k]] + facet_wrap(~variable, scales="free_x", ncol=8)
      } else if (length(unique(taxonomy_merged$variable)) >= 5) {
        pd_ggplots[[k]] <- pd_ggplots[[k]] + facet_wrap(~variable, scales="free_x", ncol=4)
      } else {
        pd_ggplots[[k]] <- pd_ggplots[[k]] + facet_wrap(~variable, scales="free_x", ncol=2)
      }
    } else {
      pd_ggplots[[k]] <- ggplot(data = taxonomy_merged, aes(value, prediction*100)) + geom_line(aes(colour=variable), size= 1) +
        scale_x_continuous(trans="log2") + labs(x="Normalized abundance", y="Prediction (% chance to be classified)", colour="Function") +
        facet_grid(class~variable, scales="free_x") + theme(legend.position="none")
    }
  }
    if ((length(unique(taxonomy_merged$variable)) >= 20 && !(modelstorun_frame$model[k] %in% continuous_variables)) | 
             (length(unique(taxonomy_merged$variable)) >= 40 && modelstorun_frame$model[k] %in% continuous_variables)) {
    ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_", modelstorun_frame$set[k], "_pd_subplot_",modelstorun_frame$model[k],".svg"), plot=pd_ggplots[[k]],
           units="mm", width=1200, height=200)
  } else if (length(unique(taxonomy_merged$variable)) >= 10) {
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_", modelstorun_frame$set[k], "_pd_subplot_",modelstorun_frame$model[k],".svg"), plot=pd_ggplots[[k]],
         units="mm", width=600, height=200)
    } else {
      ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_", modelstorun_frame$set[k], "_pd_subplot_",modelstorun_frame$model[k],".svg"), plot=pd_ggplots[[k]],
           units="mm", width=300, height=200)
    }

}


#format and save ordinations with fitted category centroids and surfaces of the gam models of environmental data
NMDS_common_Site <- data.frame(NMDS1 = ordu_common$points[,1], NMDS2 = ordu_common$points[,2], group=data.frame(sample_data(common))$Site)
ord_common_Site <- ordiellipse(ordu_common, data.frame(sample_data(common))$Site, display = "sites", 
                 kind = "se", conf = 0.95, label = T, draw = "none")
df_ell_common_Site <- data.frame()
for(g in levels(NMDS_common_Site$group)){
  df_ell_common_Site <- rbind(df_ell_common_Site, cbind(as.data.frame(with(NMDS_common_Site[NMDS_common_Site$group==g,],
                                                   vegan:::veganCovEllipse(ord_common_Site[[g]]$cov,ord_common_Site[[g]]$center,ord_common_Site[[g]]$scale))),group=g))
}

NMDS_fun_Site = data.frame(NMDS1 = ordu_func$points[,1], NMDS2 = ordu_func$points[,2], group=data.frame(sample_data(func_pruned))$Site)
ord_fun_Site <- ordiellipse(ordu_func, data.frame(sample_data(func_pruned))$Site, display = "sites", 
                            kind = "se", conf = 0.95, label = T, draw = "none")
df_ell_fun_Site <- data.frame()
for(g in levels(NMDS_fun_Site$group)){
  df_ell_fun_Site <- rbind(df_ell_fun_Site, cbind(as.data.frame(with(NMDS_fun_Site[NMDS_fun_Site$group==g,],
                                                                     vegan:::veganCovEllipse(ord_fun_Site[[g]]$cov,ord_fun_Site[[g]]$center,ord_fun_Site[[g]]$scale))),group=g))
}

if (i == 1) {
  NMDS_common_Lake = data.frame(NMDS1 = ordu_common$points[,1], NMDS2 = ordu_common$points[,2], group=data.frame(sample_data(common))$Lake)
  ord_common_Lake <- ordiellipse(ordu_common, data.frame(sample_data(common))$Lake, display = "sites", 
                                       kind = "se", conf = 0.95, label = T, draw = "none")
  df_ell_common_Lake <- data.frame()
  for(g in levels(NMDS_common_Lake$group)){
    df_ell_common_Lake <- rbind(df_ell_common_Lake, cbind(as.data.frame(with(NMDS_common_Lake[NMDS_common_Lake$group==g,],
                                                                                         vegan:::veganCovEllipse(ord_common_Lake[[g]]$cov,ord_common_Lake[[g]]$center,ord_common_Lake[[g]]$scale))),group=g))
  }


  NMDS_common_Year = data.frame(NMDS1 = ordu_common$points[,1], NMDS2 = ordu_common$points[,2], group=data.frame(sample_data(common))$Year)
  ord_common_Year <- ordiellipse(ordu_common, data.frame(sample_data(common))$Year, display = "sites", 
                                       kind = "se", conf = 0.95, label = T, draw = "none")
  df_ell_common_Year <- data.frame()
  for(g in levels(NMDS_common_Year$group)){
    df_ell_common_Year <- rbind(df_ell_common_Year, cbind(as.data.frame(with(NMDS_common_Year[NMDS_common_Year$group==g,],
                                                                                         vegan:::veganCovEllipse(ord_common_Year[[g]]$cov,ord_common_Year[[g]]$center,ord_common_Year[[g]]$scale))),group=g))
  }


ordi_common_H2S <- ordisurfs[[5]] #fetch the ordisurf object
ordi_grid_common_H2S <- ordi_common_H2S$grid #extracts the ordisurf object
ordi_a_common_H2S <- expand.grid(x = ordi_grid_common_H2S$x, y = ordi_grid_common_H2S$y) #get x and ys
ordi_a_common_H2S$z <- rescale(as.vector(ordi_grid_common_H2S$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_common_H2S <- data.frame(na.omit(ordi_a_common_H2S)) #gets rid of the nas
ordi_a_na_common_H2S$category <- "H2S"

ordi_common_Sdepth <- ordisurfs[[6]] #fetch the ordisurf object
ordi_grid_common_Sdepth <- ordi_common_Sdepth$grid #extracts the ordisurf object
ordi_a_common_Sdepth <- expand.grid(x = ordi_grid_common_Sdepth$x, y = ordi_grid_common_Sdepth$y) #get x and ys
ordi_a_common_Sdepth$z <- rescale(as.vector(ordi_grid_common_Sdepth$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_common_Sdepth <- data.frame(na.omit(ordi_a_common_Sdepth)) #gets rid of the nas
ordi_a_na_common_Sdepth$category <- "Sdepth"

ordi_common_Wdepth <- ordisurfs[[7]] #fetch the ordisurf object
ordi_grid_common_Wdepth <- ordi_common_Wdepth$grid #extracts the ordisurf object
ordi_a_common_Wdepth <- expand.grid(x = ordi_grid_common_Wdepth$x, y = ordi_grid_common_Wdepth$y) #get x and ys
ordi_a_common_Wdepth$z <- rescale(as.vector(ordi_grid_common_Wdepth$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_common_Wdepth <- data.frame(na.omit(ordi_a_common_Wdepth)) #gets rid of the nas
ordi_a_na_common_Wdepth$category <- "Wdepth"

ordi_common_pH <- ordisurfs[[8]] #fetch the ordisurf object
ordi_grid_common_pH <- ordi_common_pH$grid #extracts the ordisurf object
ordi_a_common_pH <- expand.grid(x = ordi_grid_common_pH$x, y = ordi_grid_common_pH$y) #get x and ys
ordi_a_common_pH$z <- rescale(as.vector(ordi_grid_common_pH$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_common_pH <- data.frame(na.omit(ordi_a_common_pH)) #gets rid of the nas
ordi_a_na_common_pH$category <- "pH"

ordi_common_Redox <- ordisurfs[[9]] #fetch the ordisurf object
ordi_grid_common_Redox <- ordi_common_Redox$grid #extracts the ordisurf object
ordi_a_common_Redox <- expand.grid(x = ordi_grid_common_Redox$x, y = ordi_grid_common_Redox$y) #get x and ys
ordi_a_common_Redox$z <- rescale(as.vector(ordi_grid_common_Redox$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_common_Redox <- data.frame(na.omit(ordi_a_common_Redox)) #gets rid of the nas
ordi_a_na_common_Redox$category <- "Redox"

ordi_common_O2 <- ordisurfs[[10]] #fetch the ordisurf object
ordi_grid_common_O2 <- ordi_common_O2$grid #extracts the ordisurf object
ordi_a_common_O2 <- expand.grid(x = ordi_grid_common_O2$x, y = ordi_grid_common_O2$y) #get x and ys
ordi_a_common_O2$z <- rescale(as.vector(ordi_grid_common_O2$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_common_O2 <- data.frame(na.omit(ordi_a_common_O2)) #gets rid of the nas
ordi_a_na_common_O2$category <- "O2"

ordi_common_Nseq <- ordisurfs[[11]] #fetch the ordisurf object
ordi_grid_common_Nseq <- ordi_common_Nseq$grid #extracts the ordisurf object
ordi_a_common_Nseq <- expand.grid(x = ordi_grid_common_Nseq$x, y = ordi_grid_common_Nseq$y) #get x and ys
ordi_a_common_Nseq$z <- rescale(as.vector(ordi_grid_common_Nseq$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_common_Nseq <- data.frame(na.omit(ordi_a_common_Nseq)) #gets rid of the nas
ordi_a_na_common_Nseq$category <- "Nseq"

ordisurf_NMDS_data_common <- data.frame(sample_data(common)) #there are other ways of doing this. But this is the way I do it for ease of plotting
ordisurf_NMDS_data_common$NMDS1 <- ordu_common$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
ordisurf_NMDS_data_common$NMDS2 <- ordu_common$points[ ,2]

ellipses <- c("#3f268a", #2014
              "#ffaa00", #2015
              "#910a80", #deep hole
              "#5199ff", #john's island
              "#797368", #lake hazen
              "#aac856", #skeleton lake
              "#a67e48") #snowgoose bay

contours <- list(c("#3f268a", #2014
                   "#ffaa00", #2015
                   "#910a80", #deep hole
                   #"#ff0000",
                   "#323232",
                   "#5199ff", #john's island
                   "#797368", #lake hazen
                   "#aac856", #skeleton lake
                   "#a67e48"), #snowgoose bay
                 c("#3f268a", #2014
                   "#ffaa00", #2015
                   "#910a80", #deep hole
                   "#5199ff", #john's island
                   "#797368", #lake hazen
                   #"#775500",
                   "#323232",
                   "#aac856", #skeleton lake
                   "#a67e48"), #snowgoose bay
                 c("#3f268a", #2014
                   "#ffaa00", #2015
                   "#910a80", #deep hole
                   "#5199ff", #john's island
                   "#797368", #lake hazen
                   "#aac856", #skeleton lake
                   "#a67e48", #snowgoose bay
                   #"#004cff"),
                   "#323232"),
                 c("#3f268a", #2014
                   "#ffaa00", #2015
                   "#910a80", #deep hole
                   "#5199ff", #john's island
                   "#797368", #lake hazen
                   #"#36db3c",
                   "#323232",
                   "#aac856", #skeleton lake
                   "#a67e48"), #snowgoose bay
                 c("#3f268a", #2014
                   "#ffaa00", #2015
                   "#910a80", #deep hole
                   "#5199ff", #john's island
                   "#797368", #lake hazen
                   "#323232",
                   "#aac856", #skeleton lake
                   "#a67e48"), #snowgoose bay
                 c("#3f268a", #2014
                   "#ffaa00", #2015
                   "#910a80", #deep hole
                   "#5199ff", #john's island
                   "#797368", #lake hazen
                   #"#cc44dd",
                   "#323232",
                   "#aac856", #skeleton lake
                   "#a67e48"), #snowgoose bay
                 c("#3f268a", #2014
                   "#ffaa00", #2015
                   "#910a80", #deep hole
                   "#5199ff", #john's island
                   "#797368", #lake hazen
                   "#323232",
                   "#aac856", #skeleton lake
                   "#a67e48") #snowgoose bay
                 )

plots_common <- c("ordi_a_na_common_H2S", "ordi_a_na_common_Sdepth", "ordi_a_na_common_Wdepth", 
                  "ordi_a_na_common_pH", "ordi_a_na_common_Redox", "ordi_a_na_common_O2",  "ordi_a_na_common_Nseq")


for (m in 1:length(plots_common)) {
  image <- ggplot(data = ordisurf_NMDS_data_common, aes(NMDS1, NMDS2)) + geom_text(aes(label=ordisurf_NMDS_data_common$Depth.cm, colour=Site), size= 5) +
    geom_path(data=df_ell_common_Site, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=2) +
    geom_path(data=df_ell_common_Lake, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=1) +
    geom_path(data=df_ell_common_Year, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=3) +
    geom_contour(data = eval(parse(text=plots_common[m])), aes(x = x, y = y, z = z, colour=category, alpha=..level..),bins = 20, size=1) +
    scale_color_manual(values=contours[[m]], name = "Group") + theme(legend.position = "none")
    ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_", plots_common[m], ".svg"), plot=image,
               units="mm", width=250, height=200)
}

#functionally mapped data here

NMDS_fun_Lake = data.frame(NMDS1 = ordu_func$points[,1], NMDS2 = ordu_func$points[,2], group=data.frame(sample_data(func_pruned))$Lake)
ord_fun_Lake <- ordiellipse(ordu_func, data.frame(sample_data(func_pruned))$Lake, display = "sites", 
                                     kind = "se", conf = 0.95, label = T, draw = "none")
df_ell_fun_Lake <- data.frame()
for(g in levels(NMDS_fun_Lake$group)){
  df_ell_fun_Lake <- rbind(df_ell_fun_Lake, cbind(as.data.frame(with(NMDS_fun_Lake[NMDS_fun_Lake$group==g,],
                                                                                       vegan:::veganCovEllipse(ord_fun_Lake[[g]]$cov,ord_fun_Lake[[g]]$center,ord_fun_Lake[[g]]$scale))),group=g))
}

NMDS_fun_Year = data.frame(NMDS1 = ordu_func$points[,1], NMDS2 = ordu_func$points[,2], group=data.frame(sample_data(func_pruned))$Year)
ord_fun_Year <- ordiellipse(ordu_func, data.frame(sample_data(func_pruned))$Year, display = "sites", 
                                            kind = "se", conf = 0.95, label = T, draw = "none")
df_ell_fun_Year <- data.frame()
for(g in levels(NMDS_fun_Year$group)){
  df_ell_fun_Year <- rbind(df_ell_fun_Year, cbind(as.data.frame(with(NMDS_fun_Year[NMDS_fun_Year$group==g,],
                                                                                                     vegan:::veganCovEllipse(ord_fun_Year[[g]]$cov,ord_fun_Year[[g]]$center,ord_fun_Year[[g]]$scale))),group=g))
}

#ordisurf:
ordi_fun_H2S <- ordisurfs[[16]] #fetch the ordisurf object
ordi_grid_fun_H2S <- ordi_fun_H2S$grid #extracts the ordisurf object
ordi_a_fun_H2S <- expand.grid(x = ordi_grid_fun_H2S$x, y = ordi_grid_fun_H2S$y) #get x and ys
ordi_a_fun_H2S$z <- rescale(as.vector(ordi_grid_fun_H2S$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_fun_H2S <- data.frame(na.omit(ordi_a_fun_H2S)) #gets rid of the nas
ordi_a_na_fun_H2S$category <- "H2S"

ordi_fun_Sdepth <- ordisurfs[[17]] #fetch the ordisurf object
ordi_grid_fun_Sdepth <- ordi_fun_Sdepth$grid #extracts the ordisurf object
ordi_a_fun_Sdepth <- expand.grid(x = ordi_grid_fun_Sdepth$x, y = ordi_grid_fun_Sdepth$y) #get x and ys
ordi_a_fun_Sdepth$z <- rescale(as.vector(ordi_grid_fun_Sdepth$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_fun_Sdepth <- data.frame(na.omit(ordi_a_fun_Sdepth)) #gets rid of the nas
ordi_a_na_fun_Sdepth$category <- "Sdepth"

ordi_fun_Wdepth <- ordisurfs[[18]] #fetch the ordisurf object
ordi_grid_fun_Wdepth <- ordi_fun_Wdepth$grid #extracts the ordisurf object
ordi_a_fun_Wdepth <- expand.grid(x = ordi_grid_fun_Wdepth$x, y = ordi_grid_fun_Wdepth$y) #get x and ys
ordi_a_fun_Wdepth$z <- rescale(as.vector(ordi_grid_fun_Wdepth$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_fun_Wdepth <- data.frame(na.omit(ordi_a_fun_Wdepth)) #gets rid of the nas
ordi_a_na_fun_Wdepth$category <- "Wdepth"

ordi_fun_pH <- ordisurfs[[19]] #fetch the ordisurf object
ordi_grid_fun_pH <- ordi_fun_pH$grid #extracts the ordisurf object
ordi_a_fun_pH <- expand.grid(x = ordi_grid_fun_pH$x, y = ordi_grid_fun_pH$y) #get x and ys
ordi_a_fun_pH$z <- rescale(as.vector(ordi_grid_fun_pH$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_fun_pH <- data.frame(na.omit(ordi_a_fun_pH)) #gets rid of the nas
ordi_a_na_fun_pH$category <- "pH"

ordi_fun_Redox <- ordisurfs[[20]] #fetch the ordisurf object
ordi_grid_fun_Redox <- ordi_fun_Redox$grid #extracts the ordisurf object
ordi_a_fun_Redox <- expand.grid(x = ordi_grid_fun_Redox$x, y = ordi_grid_fun_Redox$y) #get x and ys
ordi_a_fun_Redox$z <- rescale(as.vector(ordi_grid_fun_Redox$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_fun_Redox <- data.frame(na.omit(ordi_a_fun_Redox)) #gets rid of the nas
ordi_a_na_fun_Redox$category <- "Redox"

ordi_fun_O2 <- ordisurfs[[21]] #fetch the ordisurf object
ordi_grid_fun_O2 <- ordi_fun_O2$grid #extracts the ordisurf object
ordi_a_fun_O2 <- expand.grid(x = ordi_grid_fun_O2$x, y = ordi_grid_fun_O2$y) #get x and ys
ordi_a_fun_O2$z <- rescale(as.vector(ordi_grid_fun_O2$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_fun_O2 <- data.frame(na.omit(ordi_a_fun_O2)) #gets rid of the nas
ordi_a_na_fun_O2$category <- "O2"

ordi_fun_Nseq <- ordisurfs[[22]] #fetch the ordisurf object
ordi_grid_fun_Nseq <- ordi_fun_Nseq$grid #extracts the ordisurf object
ordi_a_fun_Nseq <- expand.grid(x = ordi_grid_fun_Nseq$x, y = ordi_grid_fun_Nseq$y) #get x and ys
ordi_a_fun_Nseq$z <- rescale(as.vector(ordi_grid_fun_Nseq$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_fun_Nseq <- data.frame(na.omit(ordi_a_fun_Nseq)) #gets rid of the nas
ordi_a_na_fun_Nseq$category <- "Nseq"

ordisurf_NMDS_data_fun <- data.frame(sample_data(func_pruned)) #there are other ways of doing this. But this is the way I do it for ease of plotting
ordisurf_NMDS_data_fun$NMDS1 <- ordu_func$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
ordisurf_NMDS_data_fun$NMDS2 <- ordu_func$points[ ,2]

plots_fun <- c("ordi_a_na_fun_H2S", "ordi_a_na_fun_Sdepth", "ordi_a_na_fun_Wdepth", 
                  "ordi_a_na_fun_pH", "ordi_a_na_fun_Redox", "ordi_a_na_fun_O2", "ordi_a_na_fun_Nseq")

for (m in 1:length(plots_fun)) {
  image <- ggplot(data = ordisurf_NMDS_data_fun, aes(NMDS1, NMDS2)) + geom_text(aes(label=ordisurf_NMDS_data_fun$Depth.cm, colour=Site), size= 5) +
    geom_path(data=df_ell_fun_Site, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=2) +
    geom_path(data=df_ell_fun_Lake, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=1) +
    geom_path(data=df_ell_fun_Year, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=3) +
    geom_contour(data = eval(parse(text=plots_fun[m])), aes(x = x, y = y, z = z, colour=category, alpha=..level..),bins = 20, size=1) +
    scale_color_manual(values=contours[[m]], name = "Group") + theme(legend.position = "none")
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_", plots_fun[m], ".svg"), plot=image,
         units="mm", width=250, height=200)
}

models_out <- modelstorun_frame
models_out$set <- as.vector(models_out$set)
models_out[grepl("common",models_out$set),]$set <- "OTUs"
models_out[grepl("func_pruned",models_out$set),]$set <- "FAPROTAX"
models_out$n.before <- NA
models_out$MSPE.before <- NA
models_out$CI.95.before <- NA
models_out$r2.before <- NA
models_out$Kappa.before <- NA
models_out$OOBPredErr.before <- NA
for (n in 1:nrow(models_out)) {
  models_out$n.before[n] <- models[[n]]$num.independent.variables
  if (models_out$model[n] %in% continuous_variables) {
      models_out$MSPE.before[n] <- round(all_full_model_errors[[n]][2],2)
      models_out$CI.95.before[n] <- paste0(round(all_full_model_errors[[n]][1],2),"-",round(all_full_model_errors[[n]][3],2))
      models_out$r2.before[n] <- round(models[[n]]$r.squared,3)
  } else {
      models_out$Kappa.before[n] <- round(all_full_model_errors[[n]],2)
      models_out$OOBPredErr.before[n] <- round(models[[n]]$prediction.error * 100,2)
  }
}

models_out$n.after <- NA
models_out$MSPE.after <- NA
models_out$CI.95.after <- NA
models_out$r2.after <- NA
models_out$Kappa.after <- NA
models_out$OOBPredErr.after <- NA
models_out$tax_level <- NA
for (n in 1:nrow(models_out)) {
  models_out$n.after[n] <- models_2[[n]]$num.independent.variables
  models_out$tax_level[n] <- tax_level_sets[tax_level_nums[n]]
  if (models_out$model[n] %in% continuous_variables) {
    models_out$MSPE.after[n] <- round(all_model_errors_2[[n]][2],2)
    models_out$CI.95.after[n] <- paste0(round(all_model_errors_2[[n]][1],2),"-",round(all_model_errors_2[[n]][3],2))
    models_out$r2.after[n] <- round(models_2[[n]]$r.squared,3)
  } else {
    models_out$Kappa.after[n] <- round(all_model_errors_2[[n]],2)
    models_out$OOBPredErr.after[n] <- round(models_2[[n]]$prediction.error * 100,2)
  }
}

write.csv(models_out,"D:/VirtualBox/VirtualBox Share/16S/Hazen16S/results_table.csv",quote=F)

#plot the sediment depth ordisurfs with site symbols and significant vectors for both taxonomic and functional data
physchem_vectors_common <- physchem_vectors[physchem_vectors$set %in% "OTUs",]
plot(ordu_common)
multiplier_common <- ordiArrowMul(ef_common)
image <- ggplot(data = ordisurf_NMDS_data_common, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(shape=Lake, colour=Site), size= 5) + coord_fixed() +
  geom_path(data=df_ell_common_Site, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=2) +
  geom_path(data=df_ell_common_Lake, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=1) +
  geom_path(data=df_ell_common_Year, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=3) +
  geom_contour(data = eval(parse(text=plots_common[2])), aes(x = x, y = y, z = z, colour=category, alpha=..level..),bins = 10, size=1) +
  geom_segment(data=physchem_vectors_common, aes(x=0, xend=NMDS1*multiplier_common, y=0, yend=NMDS2*multiplier_common), arrow = arrow(length = unit(0.5, "cm")), colour="black", size=1) +
  geom_text(data=physchem_vectors_common, aes(x=NMDS1*multiplier_common+multiplier_common*0.1, y=NMDS2*multiplier_common+multiplier_common*0.1, label=model),size=5) + 
  annotate("text", label = paste0("Stress = ", round(ordu_common$stress,3)), x = min(ordisurf_NMDS_data_common$NMDS1)*0.8, y = min(ordisurf_NMDS_data_common$NMDS2)*0.9, size = 6) + 
  scale_color_manual(values=contours[[2]], name = "Group") + theme(legend.position = "none")
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_ordination_vectors_tax.svg"), plot=image,
       units="mm", width=300, height=300)

physchem_vectors_func <- physchem_vectors[physchem_vectors$set %in% "FAPROTAX",]
plot(ordu_func)
multiplier_func <- ordiArrowMul(ef_func)
image <- ggplot(data = ordisurf_NMDS_data_fun, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(shape=Lake, colour=Site), size= 5) + coord_fixed() +
  geom_path(data=df_ell_fun_Site, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=2) +
  geom_path(data=df_ell_fun_Lake, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=1) +
  geom_path(data=df_ell_fun_Year, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=3) +
  geom_contour(data = eval(parse(text=plots_fun[2])), aes(x = x, y = y, z = z, colour=category, alpha=..level..),bins = 10, size=1) +
  geom_segment(data=physchem_vectors_func, aes(x=0, xend=NMDS1*multiplier_func, y=0, yend=NMDS2*multiplier_func), arrow = arrow(length = unit(0.5, "cm")), colour="black", size=1) +
  geom_text(data=physchem_vectors_func, aes(x=NMDS1*multiplier_func+multiplier_func*0.1, y=NMDS2*multiplier_func+multiplier_func*0.1, label=model),size=5) +
  annotate("text", label = paste0("Stress = ", round(ordu_func$stress,3)), x = min(ordisurf_NMDS_data_fun$NMDS1)*0.8, y = min(ordisurf_NMDS_data_fun$NMDS2)*0.9, size = 6) + 
  scale_color_manual(values=contours[[2]], name = "Group") + theme(legend.position = "none")
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_ordination_vectors_fun.svg"), plot=image,
       units="mm", width=300, height=300)


#and run the Summer 2015 data sets

} else {
  ordi_common_Sdepth <- ordisurfs[[3]] #fetch the ordisurf object
  ordi_grid_common_Sdepth <- ordi_common_Sdepth$grid #extracts the ordisurf object
  ordi_a_common_Sdepth <- expand.grid(x = ordi_grid_common_Sdepth$x, y = ordi_grid_common_Sdepth$y) #get x and ys
  ordi_a_common_Sdepth$z <- rescale(as.vector(ordi_grid_common_Sdepth$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_common_Sdepth <- data.frame(na.omit(ordi_a_common_Sdepth)) #gets rid of the nas
  ordi_a_na_common_Sdepth$category <- "Sdepth"
  
  ordi_common_pH <- ordisurfs[[4]] #fetch the ordisurf object
  ordi_grid_common_pH <- ordi_common_pH$grid #extracts the ordisurf object
  ordi_a_common_pH <- expand.grid(x = ordi_grid_common_pH$x, y = ordi_grid_common_pH$y) #get x and ys
  ordi_a_common_pH$z <- rescale(as.vector(ordi_grid_common_pH$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_common_pH <- data.frame(na.omit(ordi_a_common_pH)) #gets rid of the nas
  ordi_a_na_common_pH$category <- "pH"
  
  ordi_common_O2 <- ordisurfs[[5]] #fetch the ordisurf object
  ordi_grid_common_O2 <- ordi_common_O2$grid #extracts the ordisurf object
  ordi_a_common_O2 <- expand.grid(x = ordi_grid_common_O2$x, y = ordi_grid_common_O2$y) #get x and ys
  ordi_a_common_O2$z <- rescale(as.vector(ordi_grid_common_O2$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_common_O2 <- data.frame(na.omit(ordi_a_common_O2)) #gets rid of the nas
  ordi_a_na_common_O2$category <- "O2"
  
  ordi_common_NO3 <- ordisurfs[[6]] #fetch the ordisurf object
  ordi_grid_common_NO3 <- ordi_common_NO3$grid #extracts the ordisurf object
  ordi_a_common_NO3 <- expand.grid(x = ordi_grid_common_NO3$x, y = ordi_grid_common_NO3$y) #get x and ys
  ordi_a_common_NO3$z <- rescale(as.vector(ordi_grid_common_NO3$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_common_NO3 <- data.frame(na.omit(ordi_a_common_NO3)) #gets rid of the nas
  ordi_a_na_common_NO3$category <- "NO3"
  
  ordi_common_Cl <- ordisurfs[[7]] #fetch the ordisurf object
  ordi_grid_common_Cl <- ordi_common_Cl$grid #extracts the ordisurf object
  ordi_a_common_Cl <- expand.grid(x = ordi_grid_common_Cl$x, y = ordi_grid_common_Cl$y) #get x and ys
  ordi_a_common_Cl$z <- rescale(as.vector(ordi_grid_common_Cl$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_common_Cl <- data.frame(na.omit(ordi_a_common_Cl)) #gets rid of the nas
  ordi_a_na_common_Cl$category <- "Chloride"
  
  ordi_common_SO42 <- ordisurfs[[8]] #fetch the ordisurf object
  ordi_grid_common_SO42 <- ordi_common_SO42$grid #extracts the ordisurf object
  ordi_a_common_SO42 <- expand.grid(x = ordi_grid_common_SO42$x, y = ordi_grid_common_SO42$y) #get x and ys
  ordi_a_common_SO42$z <- rescale(as.vector(ordi_grid_common_SO42$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_common_SO42 <- data.frame(na.omit(ordi_a_common_SO42)) #gets rid of the nas
  ordi_a_na_common_SO42$category <- "SO42"
  
  ordi_common_Nseq <- ordisurfs[[9]] #fetch the ordisurf object
  ordi_grid_common_Nseq <- ordi_common_Nseq$grid #extracts the ordisurf object
  ordi_a_common_Nseq <- expand.grid(x = ordi_grid_common_Nseq$x, y = ordi_grid_common_Nseq$y) #get x and ys
  ordi_a_common_Nseq$z <- rescale(as.vector(ordi_grid_common_Nseq$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_common_Nseq <- data.frame(na.omit(ordi_a_common_Nseq)) #gets rid of the nas
  ordi_a_na_common_Nseq$category <- "Nseq"
  
  ordisurf_NMDS_data_common <- data.frame(sample_data(common)) #there are other ways of doing this. But this is the way I do it for ease of plotting
  ordisurf_NMDS_data_common$NMDS1 <- ordu_common$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
  ordisurf_NMDS_data_common$NMDS2 <- ordu_common$points[ ,2]

  ellipses <- c("#229373", #pond1
                "#aac856") #skeleton lake
  
  contours <- list(c("#229373", #pond1
                     #"#775500",
                     "#323232",
                     "#aac856"), #skeleton lake
                   c(#"#cc44dd",
                     "#323232",
                     "#229373", #pond1
                     "#aac856"), #skeleton lake
                   c(#"#36db3c",
                     "#323232",
                     "#229373", #pond1
                     "#aac856"), #skeleton lake
                   c(#"#e0b702",
                     "#323232",
                     "#229373", #pond1
                     "#aac856"), #skeleton lake
                   c(#"#60ffbf",
                     "#323232",
                     "#229373", #pond1
                     "#aac856"), #skeleton lake
                   c("#229373", #pond1
                     "#aac856", #skeleton lake
                     #"#99163b")
                     "#323232"),
                   c("#323232",
                     "#229373", #pond1
                     "#aac856") #skeleton lake
  )
  
  plots_common <- c("ordi_a_na_common_Sdepth", "ordi_a_na_common_O2", "ordi_a_na_common_pH",
                    "ordi_a_na_common_NO3", "ordi_a_na_common_Cl", "ordi_a_na_common_SO42", "ordi_a_na_common_Nseq")
  
  
  for (m in 1:length(plots_common)) {
    image <- ggplot(data = ordisurf_NMDS_data_common, aes(NMDS1, NMDS2)) + geom_text(aes(label=ordisurf_NMDS_data_common$Depth.cm, colour=Site), size= 5) +
      geom_path(data=df_ell_common_Site, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=2) +
      geom_contour(data = eval(parse(text=plots_common[m])), aes(x = x, y = y, z = z, colour=category, alpha=..level..),bins = 20, size=1) +
      scale_color_manual(values=contours[[m]], name = "Group") + theme(legend.position = "none")
    ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_", plots_common[m], ".svg"), plot=image,
           units="mm", width=250, height=200)
  }
  
  #functionally mapped data here
  
  ordi_fun_Sdepth <- ordisurfs[[12]] #fetch the ordisurf object
  ordi_grid_fun_Sdepth <- ordi_fun_Sdepth$grid #extracts the ordisurf object
  ordi_a_fun_Sdepth <- expand.grid(x = ordi_grid_fun_Sdepth$x, y = ordi_grid_fun_Sdepth$y) #get x and ys
  ordi_a_fun_Sdepth$z <- rescale(as.vector(ordi_grid_fun_Sdepth$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_fun_Sdepth <- data.frame(na.omit(ordi_a_fun_Sdepth)) #gets rid of the nas
  ordi_a_na_fun_Sdepth$category <- "Sdepth"
  
  ordi_fun_pH <- ordisurfs[[13]] #fetch the ordisurf object
  ordi_grid_fun_pH <- ordi_fun_pH$grid #extracts the ordisurf object
  ordi_a_fun_pH <- expand.grid(x = ordi_grid_fun_pH$x, y = ordi_grid_fun_pH$y) #get x and ys
  ordi_a_fun_pH$z <- rescale(as.vector(ordi_grid_fun_pH$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_fun_pH <- data.frame(na.omit(ordi_a_fun_pH)) #gets rid of the nas
  ordi_a_na_fun_pH$category <- "pH"
  
  ordi_fun_O2 <- ordisurfs[[14]] #fetch the ordisurf object
  ordi_grid_fun_O2 <- ordi_fun_O2$grid #extracts the ordisurf object
  ordi_a_fun_O2 <- expand.grid(x = ordi_grid_fun_O2$x, y = ordi_grid_fun_O2$y) #get x and ys
  ordi_a_fun_O2$z <- rescale(as.vector(ordi_grid_fun_O2$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_fun_O2 <- data.frame(na.omit(ordi_a_fun_O2)) #gets rid of the nas
  ordi_a_na_fun_O2$category <- "O2"
  
  ordi_fun_NO3 <- ordisurfs[[15]] #fetch the ordisurf object
  ordi_grid_fun_NO3 <- ordi_fun_NO3$grid #extracts the ordisurf object
  ordi_a_fun_NO3 <- expand.grid(x = ordi_grid_fun_NO3$x, y = ordi_grid_fun_NO3$y) #get x and ys
  ordi_a_fun_NO3$z <- rescale(as.vector(ordi_grid_fun_NO3$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_fun_NO3 <- data.frame(na.omit(ordi_a_fun_NO3)) #gets rid of the nas
  ordi_a_na_fun_NO3$category <- "NO3"
  
  ordi_fun_Cl <- ordisurfs[[16]] #fetch the ordisurf object
  ordi_grid_fun_Cl <- ordi_fun_Cl$grid #extracts the ordisurf object
  ordi_a_fun_Cl <- expand.grid(x = ordi_grid_fun_Cl$x, y = ordi_grid_fun_Cl$y) #get x and ys
  ordi_a_fun_Cl$z <- rescale(as.vector(ordi_grid_fun_Cl$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_fun_Cl <- data.frame(na.omit(ordi_a_fun_Cl)) #gets rid of the nas
  ordi_a_na_fun_Cl$category <- "Chloride"
  
  ordi_fun_SO42 <- ordisurfs[[17]] #fetch the ordisurf object
  ordi_grid_fun_SO42 <- ordi_fun_SO42$grid #extracts the ordisurf object
  ordi_a_fun_SO42 <- expand.grid(x = ordi_grid_fun_SO42$x, y = ordi_grid_fun_SO42$y) #get x and ys
  ordi_a_fun_SO42$z <- rescale(as.vector(ordi_grid_fun_SO42$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_fun_SO42 <- data.frame(na.omit(ordi_a_fun_SO42)) #gets rid of the nas
  ordi_a_na_fun_SO42$category <- "SO42"
  
  ordi_fun_Nseq <- ordisurfs[[18]] #fetch the ordisurf object
  ordi_grid_fun_Nseq <- ordi_fun_Nseq$grid #extracts the ordisurf object
  ordi_a_fun_Nseq <- expand.grid(x = ordi_grid_fun_Nseq$x, y = ordi_grid_fun_Nseq$y) #get x and ys
  ordi_a_fun_Nseq$z <- rescale(as.vector(ordi_grid_fun_Nseq$z), to=c(0,100)) #unravel the matrix for the z scores
  ordi_a_na_fun_Nseq <- data.frame(na.omit(ordi_a_fun_Nseq)) #gets rid of the nas
  ordi_a_na_fun_Nseq$category <- "Nseq"

  
  ordisurf_NMDS_data_fun <- data.frame(sample_data(func_pruned)) #there are other ways of doing this. But this is the way I do it for ease of plotting
  ordisurf_NMDS_data_fun$NMDS1 <- ordu_func$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
  ordisurf_NMDS_data_fun$NMDS2 <- ordu_func$points[ ,2]
  
  plots_fun <- c("ordi_a_na_fun_Sdepth", "ordi_a_na_fun_O2", "ordi_a_na_fun_pH",
                 "ordi_a_na_fun_NO3", "ordi_a_na_fun_Cl", "ordi_a_na_fun_SO42", "ordi_a_na_fun_Nseq")
  
  for (m in 1:length(plots_common)) {
    image <- ggplot(data = ordisurf_NMDS_data_fun, aes(NMDS1, NMDS2)) + geom_text(aes(label=ordisurf_NMDS_data_fun$Depth.cm, colour=Site), size= 5) +
      geom_path(data=df_ell_fun_Site, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=2) +
      geom_contour(data = eval(parse(text=plots_fun[m])), aes(x = x, y = y, z = z, colour=category, alpha=..level..),bins = 20, size=1) +
      scale_color_manual(values=contours[[m]], name = "Group") + theme(legend.position = "none")
    ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_", plots_fun[m], ".svg"), plot=image,
           units="mm", width=250, height=200)
  }
  
  models_out <- modelstorun_frame
  models_out$set <- as.vector(models_out$set)
  models_out[grepl("common",models_out$set),]$set <- "OTUs"
  models_out[grepl("func_pruned",models_out$set),]$set <- "FAPROTAX"
  models_out$n.before <- NA
  models_out$MSPE.before <- NA
  models_out$CI.95.before <- NA
  models_out$r2.before <- NA
  models_out$Kappa.before <- NA
  models_out$OOBPredErr.before <- NA
  for (n in 1:nrow(models_out)) {
    models_out$n.before[n] <- models[[n]]$num.independent.variables
    if (models_out$model[n] %in% continuous_variables) {
      models_out$MSPE.before[n] <- round(all_full_model_errors[[n]][2],2)
      models_out$CI.95.before[n] <- paste0(round(all_full_model_errors[[n]][1],2),"-",round(all_full_model_errors[[n]][3],2))
      models_out$r2.before[n] <- round(models[[n]]$r.squared,3)
    } else {
      models_out$Kappa.before[n] <- round(all_full_model_errors[[n]],2)
      models_out$OOBPredErr.before[n] <- round(models[[n]]$prediction.error * 100,2)
    }
  }
  
  models_out$n.after <- NA
  models_out$MSPE.after <- NA
  models_out$CI.95.after <- NA
  models_out$r2.after <- NA
  models_out$Kappa.after <- NA
  models_out$OOBPredErr.after <- NA
  models_out$tax_level <- NA
  for (n in 1:nrow(models_out)) {
    models_out$n.after[n] <- models_2[[n]]$num.independent.variables
    models_out$tax_level[n] <- tax_level_sets[tax_level_nums[n]]
    if (models_out$model[n] %in% continuous_variables) {
      models_out$MSPE.after[n] <- round(all_model_errors_2[[n]][2],2)
      models_out$CI.95.after[n] <- paste0(round(all_model_errors_2[[n]][1],2),"-",round(all_model_errors_2[[n]][3],2))
      models_out$r2.after[n] <- round(models_2[[n]]$r.squared,3)
    } else {
      models_out$Kappa.after[n] <- round(all_model_errors_2[[n]],2)
      models_out$OOBPredErr.after[n] <- round(models_2[[n]]$prediction.error * 100,2)
    }
  }
  write.csv(models_out,paste0("D:/VirtualBox/VirtualBox Share/16S/hazensummer/", dataset_names[i], "_results_table.csv"),quote=F)
  
  #plot the sediment depth ordisurfs with site symbols and significant vectors for both taxonomic and functional data
  physchem_vectors_common <- physchem_vectors[physchem_vectors$set %in% "OTUs",]
  plot(ordu_common)
  multiplier_common <- ordiArrowMul(ef_common)
  image <- ggplot(data = ordisurf_NMDS_data_common, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(colour=Site), size= 5) + coord_fixed(ratio=1, expand=T) +
    geom_path(data=df_ell_common_Site, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=2) +
    geom_contour(data = eval(parse(text=plots_common[1])), aes(x = x, y = y, z = z, colour=category, alpha=..level..),bins = 10, size=1) +
    geom_segment(data=physchem_vectors_common, aes(x=0, xend=NMDS1*multiplier_common, y=0, yend=NMDS2*multiplier_common), arrow = arrow(length = unit(0.5, "cm")), colour="black", size=1) +
    geom_text(data=physchem_vectors_common, aes(x=NMDS1*multiplier_common+multiplier_common*0.1, y=NMDS2*multiplier_common+multiplier_common*0.1, label=model),size=5) +
    annotate("text", label = paste0("Stress = ", round(ordu_common$stress,3)), x = min(ordisurf_NMDS_data_common$NMDS1)*0.8, y = min(ordisurf_NMDS_data_common$NMDS2)*0.9, size = 6) +
    scale_color_manual(values=contours[[1]], name = "Group") + theme(legend.position = "none") 
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_ordination_vectors_tax.svg"), plot=image,
         units="mm", width=250, height=200)
  
  physchem_vectors_func <- physchem_vectors[physchem_vectors$set %in% "FAPROTAX",]
  plot(ordu_func)
  multiplier_func <- ordiArrowMul(ef_func)
  image <- ggplot(data = ordisurf_NMDS_data_fun, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(colour=Site), size= 5) + coord_fixed(ratio=1, expand=T) +
    geom_path(data=df_ell_fun_Site, aes(x=NMDS1, y=NMDS2, colour=group), size=1.25, linetype=2) +
    geom_contour(data = eval(parse(text=plots_fun[1])), aes(x = x, y = y, z = z, colour=category, alpha=..level..),bins = 10, size=1) +
    geom_segment(data=physchem_vectors_func, aes(x=0, xend=NMDS1*multiplier_func, y=0, yend=NMDS2*multiplier_func), arrow = arrow(length = unit(0.5, "cm")), colour="black", size=1) +
    geom_text(data=physchem_vectors_func, aes(x=NMDS1*multiplier_func+multiplier_func*0.1, y=NMDS2*multiplier_func+multiplier_func*0.1, label=model),size=5) + 
    annotate("text", label = paste0("Stress = ", round(ordu_func$stress,3)), x = min(ordisurf_NMDS_data_fun$NMDS1)*0.8, y = min(ordisurf_NMDS_data_fun$NMDS2)*0.9, size = 6) +
    scale_color_manual(values=contours[[1]], name = "Group") + theme(legend.position = "none")
  ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/16S/Plots/", dataset_names[i], "_ordination_vectors_fun.svg"), plot=image,
         units="mm", width=250, height=200)
 
}
#garbage collection to free up RAM
rm(comp_classify)
gc()
}

#write out the richness and diversity model data, and mantel test summaries
richness_list <- do.call("rbind", richness_list)
write.csv(richness_list,"D:/VirtualBox/VirtualBox Share/16S/richness_list.csv",quote=F)
diversity_list <- do.call("rbind", diversity_list)
write.csv(diversity_list,"D:/VirtualBox/VirtualBox Share/16S/diversity_list.csv",quote=F)
mantel_list <- do.call("rbind", mantels)
write.csv(mantel_list,"D:/VirtualBox/VirtualBox Share/16S/mantel_list.csv",quote=F)

#write out high level description of phyla abundances and metadata
for (i in 1:length(phyloseqs)) {
unknown_OTU <- rownames(tax_table(phyloseqs[[i]])[tax_table(phyloseqs[[i]])[,2] %in% NA,2])
phyla <- merge_taxa(phyloseqs[[i]], unknown_OTU)
phyla_prot <- phyla
if(any(tax_table(phyla)[,2] %in% "Proteobacteria")) {
  replacement_proteobacteria <- tax_table(phyla)[tax_table(phyla)[,2] %in% "Proteobacteria",3]
  replacement_proteobacteria[replacement_proteobacteria %in% NA] <- "Proteobacteria"
  tax_table(phyla)[tax_table(phyla)[,2] %in% "Proteobacteria",2] <- replacement_proteobacteria
}
phyla <- tax_glom(phyla, taxrank = "Phylum")
phyla_prot <- tax_glom(phyla_prot, taxrank = "Phylum")
least_abundant_phyla <- names(taxa_sums(phyla))[taxa_sums(phyla)/sum(taxa_sums(phyla))*100 < 1]
phyla <- merge_taxa(phyla, least_abundant_phyla, 1)
tax_table(phyla)[tax_table(phyla)[,2] %in% NA,2] <- "Other / Unknown"
tax_table(phyla_prot)[tax_table(phyla_prot)[,2] %in% NA,2] <- "Other / Unknown"

phyla <- transform_sample_counts(phyla, function(x) 100 * x/sum(x))
phyla_prot <- transform_sample_counts(phyla_prot, function(x) 100 * x/sum(x))
otu_table_phyla <- data.frame(t(otu_table(phyla)))
otu_table_phyla_prot <- data.frame(t(otu_table(phyla_prot)))
colnames(otu_table_phyla) <- c(tax_table(phyla)[,2])
colnames(otu_table_phyla_prot) <- c(tax_table(phyla_prot)[,2])
phyla_levels <- phyla_levels_all[phyla_levels_all %in% colnames(otu_table_phyla)]
phyla_levels_prot <- phyla_levels_all[phyla_levels_all %in% colnames(otu_table_phyla_prot)]

otu_table_phyla <- melt(as.matrix(otu_table_phyla), varnames=c("Sample", "Phylum"), value.name="Abundance")
otu_table_phyla_prot <- melt(as.matrix(otu_table_phyla_prot), varnames=c("Sample", "Phylum"), value.name="Abundance")
otu_table_phyla$Phylum <- factor(otu_table_phyla$Phylum, levels = phyla_levels)
otu_table_phyla_prot$Phylum <- factor(otu_table_phyla_prot$Phylum, levels = phyla_levels_prot)
#means
mean_frame <- aggregate(Abundance ~ Phylum, otu_table_phyla, mean)
mean_frame_prot <- aggregate(Abundance ~ Phylum, otu_table_phyla_prot, mean)
se_frame <- aggregate(Abundance ~ Phylum, otu_table_phyla, se)
se_frame_prot <- aggregate(Abundance ~ Phylum, otu_table_phyla_prot, se)
result_frame <- merge(mean_frame,se_frame,by="Phylum")
result_frame_prot <- merge(mean_frame_prot,se_frame_prot,by="Phylum")
colnames(result_frame) <- c("Phylum","mean","se")
colnames(result_frame_prot) <- c("Phylum","mean","se")
result_frame <- result_frame[order(result_frame$mean),]
result_frame_prot <- result_frame_prot[order(result_frame_prot$mean),]
result_frame$Lake <- NA 
result_frame_prot$Lake <- NA
lakes_frame <- otu_table_phyla
if (i == 1) {lakes_frame$Lake <- "Lake Hazen"} else {lakes_frame$Lake <- "Pond1"}
lakes_frame$Lake[grep("Skeleton", lakes_frame$Sample)] <- "Skeleton Lake"
lakes_frame_prot <- otu_table_phyla_prot
if (i == 1) {lakes_frame_prot$Lake <- "Lake Hazen"} else {lakes_frame_prot$Lake <- "Pond1"}
lakes_frame_prot$Lake[grep("Skeleton", lakes_frame_prot$Sample)] <- "Skeleton Lake"
mean_frame_lake <- aggregate(Abundance ~ Phylum + Lake, lakes_frame, mean)
mean_frame_lake_prot <- aggregate(Abundance ~ Phylum + Lake, lakes_frame_prot, mean)
se_frame_lake <- aggregate(Abundance ~ Phylum + Lake, lakes_frame, se)
se_frame_lake_prot <- aggregate(Abundance ~ Phylum + Lake, lakes_frame_prot, se)
result_frame_lake <- cbind(mean_frame_lake[-3], mean=mean_frame_lake$Abundance, se=se_frame_lake$Abundance)
result_frame_lake_prot <- cbind(mean_frame_lake_prot[-3], mean=mean_frame_lake_prot$Abundance, se=se_frame_lake_prot$Abundance)
result_frame_lake <- result_frame_lake[order(result_frame_lake$Lake,result_frame_lake$mean),]
result_frame_lake_prot <- result_frame_lake_prot[order(result_frame_lake_prot$Lake,result_frame_lake_prot$mean),]
result_frame <- rbind(result_frame,result_frame_lake)
result_frame_prot <- rbind(result_frame_prot,result_frame_lake_prot)

write.csv(result_frame, paste0("D:/VirtualBox/VirtualBox Share/16S/", dataset_names[i],"_relative_abundances.csv"),quote=F)
write.csv(result_frame_prot, paste0("D:/VirtualBox/VirtualBox Share/16S/", dataset_names[i],"_relative_abundances_prot.csv"),quote=F)
#metadata
write.csv(sample_data(phyloseqs[[i]]), paste0("D:/VirtualBox/VirtualBox Share/16S/", dataset_names[i],"_metadata.csv"),quote=F)
}


# #garbage collection
# save.image("C:/Users/Matti/Documents/Hazen.RData")
# rm()
# gc()
# load("C:/Users/Matti/Documents/Hazen.RData")

#post hoc, check how analyzing pond1 and skeleton samples separately affects mantel tests
posthoc_mantels <- list(NULL)
for (i in 2:3) {
if (i==2){index <- c(1:4)} else {index <- c(5:8)}
common_skel_only <- subset_samples(phyloseqs[[i]], Site %in% "Skeleton Lake")
common_skel_only <- prune_taxa(taxa_sums(common_skel_only) > 0, common_skel_only)
common_skel_distance <- DPCoA(common_skel_only)
common_skel_sample_data <- data.frame(sample_data(common_skel_only)[,which(colnames(sample_data(common_skel_only)) %in% continuous_variables)])
common_skel_env_distance <- vegdist(common_skel_sample_data, method="euclidean")
common_skel_mantel <- mantel(common_skel_distance$RaoDis, common_skel_env_distance, method="pearson", permutations = 10000)
posthoc_mantels[[index[1]]] <- common_skel_mantel

fun_skel_only <- subset_samples(func_data[[i]], Site %in% "Skeleton Lake")
fun_skel_only <- prune_taxa(taxa_sums(fun_skel_only) > 0, fun_skel_only)
fun_skel_only_matrix <- wisconsin(sqrt(veganifyOTU(fun_skel_only)))
fun_skel_distance <- vegdist(fun_skel_only_matrix, distance = "bray")
fun_skel_sample_data <- data.frame(sample_data(fun_skel_only)[,which(colnames(sample_data(fun_skel_only)) %in% continuous_variables)])
fun_skel_env_distance <- vegdist(fun_skel_sample_data, method="euclidean")
fun_skel_mantel <- mantel(fun_skel_distance, fun_skel_env_distance, method="pearson", permutations = 10000)
posthoc_mantels[[index[2]]] <- fun_skel_mantel

common_pond1_only <- subset_samples(phyloseqs[[i]], Site %in% "Pond1")
common_pond1_only <- prune_taxa(taxa_sums(common_pond1_only) > 0, common_pond1_only)
common_pond1_distance <- DPCoA(common_pond1_only)
common_pond1_sample_data <- data.frame(sample_data(common_pond1_only)[,which(colnames(sample_data(common_pond1_only)) %in% continuous_variables)])
common_pond1_env_distance <- vegdist(common_pond1_sample_data, method="euclidean")
common_pond1_mantel <- mantel(common_pond1_distance$RaoDis, common_pond1_env_distance, method="pearson", permutations = 10000)
posthoc_mantels[[index[3]]] <- common_pond1_mantel

fun_pond1_only <- subset_samples(func_data[[i]], Site %in% "Pond1")
fun_pond1_only <- prune_taxa(taxa_sums(fun_pond1_only) > 0, fun_pond1_only)
fun_pond1_only_matrix <- wisconsin(sqrt(veganifyOTU(fun_pond1_only)))
fun_pond1_distance <- vegdist(fun_pond1_only_matrix, distance = "bray")
fun_pond1_sample_data <- data.frame(sample_data(fun_pond1_only)[,which(colnames(sample_data(fun_pond1_only)) %in% continuous_variables)])
fun_pond1_env_distance <- vegdist(fun_pond1_sample_data, method="euclidean")
fun_pond1_mantel <- mantel(fun_pond1_distance, fun_pond1_env_distance, method="pearson", permutations = 10000)
posthoc_mantels[[index[4]]] <- fun_pond1_mantel
}
posthoc_data <- data.frame(Lake = c(rep(c(rep("Skeleton Lake", 2), rep("Pond1", 2)),2)),
                           Set = c(rep(c("common", "fun", "common", "fun"), 2)),
                           R2 = sapply(posthoc_mantels, "[[", 3),
                           bonf.pvalues = p.adjust(c(sapply(posthoc_mantels, "[[", 4)), method = "bonferroni"))

#which taxa were assigned to the functional group of sulfur_respiration in spring 2014/2015?

assoc <- read.csv("D:/VirtualBox/VirtualBox Share/16S/hazen_OTU_func_association.csv", sep = "\t")
sulf <- assoc[c(1,which(colnames(assoc) %in% "sulfur_respiration"))]
sulf <- sulf[which(rowSums(sulf[2]) > 0),]
sulf_respirers <- subset_taxa(phyloseqs[[1]], taxa_names(phyloseqs[[1]]) %in% sulf$record)
plot_bar(sulf_respirers, fill = "Family")
plot_bar(sulf_respirers, fill = "Genus")

#which taxa were assigned to the functional group of mercury_methylation in summer 2015 archaea?

assoc2 <- read.csv("D:/VirtualBox/VirtualBox Share/16S/hazensummer_a_OTU_func_association.csv", sep = "\t")
merc <- assoc2[c(1,which(colnames(assoc2) %in% "mercury_methylation"))]
merc <- merc[which(rowSums(merc[2]) > 0),]
merc_methylators <- subset_taxa(phyloseqs[[2]], taxa_names(phyloseqs[[2]]) %in% merc$record)
plot_bar(merc_methylators, fill = "Species")

assoc3 <- read.csv("D:/VirtualBox/VirtualBox Share/16S/hazensummer_b_OTU_func_association.csv", sep = "\t")


#extra tables for reviews
if (identical(rownames(data.frame(tax_table(phyloseqs[[1]]))), as.vector(assoc$record))) {
  assoc <- cbind(data.frame(tax_table(phyloseqs[[1]])), assoc[-1])
}

arm_hazen <- assoc[,c(1:7,grep("intracellular_parasites", colnames(assoc)))]
arm_hazen <- arm_hazen[which(arm_hazen[,ncol(arm_hazen)] != 0),]

if (identical(rownames(data.frame(tax_table(phyloseqs[[2]]))), as.vector(assoc2$record))) {
  assoc2 <- cbind(data.frame(tax_table(phyloseqs[[2]])), assoc2[-1])
}
methanogens_hazensummer_a <- assoc2[,c(1:7,grep("methanogenesis", colnames(assoc2)))]
methanogens_hazensummer_a <- methanogens_hazensummer_a[which(rowSums(methanogens_hazensummer_a[,8:ncol(methanogens_hazensummer_a)]) != 0),]

methanogens_hazensummer_a <- assoc2[,c(1:7,grep("methanogenesis", colnames(assoc2)))]
methanogens_hazensummer_a <- methanogens_hazensummer_a[which(rowSums(methanogens_hazensummer_a[,8:ncol(methanogens_hazensummer_a)]) != 0),]

if (identical(rownames(data.frame(tax_table(phyloseqs[[3]]))), as.vector(assoc3$record))) {
  assoc3 <- cbind(data.frame(tax_table(phyloseqs[[3]])), assoc3[-1])
}
methanotrophs_hazensummer_b <- assoc3[,c(1:7,grep("methanotrophy", colnames(assoc3)))]
methanotrophs_hazensummer_b <- methanotrophs_hazensummer_b[which(methanotrophs_hazensummer_b[,ncol(methanotrophs_hazensummer_b)] != 0),]

#save OTU tables with taxonomy
for (i in 1:3) {
  if (i > 1) {
    otu_table <- data.frame(otu_table(phyloseqs[[i]]))
    colnames(otu_table) <- sample_names(phyloseqs[[i]])
    setcolorder(otu_table, sample_levels_summer[-9])
    } else {
      otu_table <- otu_table(phyloseqs[[i]])
      }
  write.csv(cbind(tax_table(phyloseqs[[i]]), otu_table), paste0("D:/VirtualBox/VirtualBox Share/16S/", dataset_names[i],"_otu_tax_table.csv"),quote=F)
  }
