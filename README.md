## Contents

This repository contains supplementary code, databases and scripts for the study "Physicochemical Drivers of Microbial Community Structure in Sediments of Lake Hazen, Nunavut, Canada" by Matti O. Ruuskanen, Kyra A. St.Pierre, Vincent L. St.Louis, Stéphane Aris-Brosou & Alexandre J. Poulain (https://doi.org/10.3389/fmicb.2018.01138).

## Included files

1. Sequence handling shell scripts
	- Scripts that were run on the CAC cluster to produce the primary data from Illumina read files
	- One shell script for each of the three data sets (Spring 2014/2015, Summer 2015 Archaea, Summer 2015 Bacteria)
	- Format: *Dataset*_Sequence_handling_shell_script.txt

2. Accessory R scripts (used by the sequence handling scripts)
	- swarm_construct_otu_table.R
		* Constructs an OTU table from output of Swarm 2.0
	- uc_to_OTU_table.R
		* Constructs on OTU table from .uc cluster files
		
3. Data analysis scripts
	- Data_analysis_scripts.R
	- Scripts that use the primary data to produce all the figures, tables, and other results in the manuscript
	- Also includes code run on CAC for FAPROTAX functional mapping (as comments)

4. Custom FAPROTAX databases
	- The database that was used to map OTU abundances into functions by taxonomy
	- All the changes made to the original database can be found as comments at the beginning of the database
5. Mapping files for demultiplexing
	- The shell scripts require four mapping files which include information on primers for demultiplexing of the reads.
	- qhazen_map_meta.txt
		* 2014 Spring samples from Lake Hazen
	- rhazen_map_meta.txt
		* 2015 Spring samples from Lake Hazen and Skeleton Lake
	- hazensummer_a_map_meta.txt
		* 2015 Summer samples from Pond1 and Skeleton Lake, archaeal reads
	- hazensummer_b_map_meta.txt
		* 2015 Summer samples from Pond1 and Skeleton Lake, bacterial reads
