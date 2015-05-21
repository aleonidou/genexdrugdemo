# ============================================ #
# ICR BIC GDSC workshop
# https://github.com/DrJCampbell/genexdrugdemo
# ============================================ #

# change the next line to the directory where your files are
setwd("./")

# file with gene mutations and copy number calls
mutation_data <- read.table(
	file="gdsc_CNV_mutation_predictors.txt",
	sep="\t",
	header=TRUE,
	row.names=1
	)

# drug IC50 values
drug_data <- read.table(
	file="gdsc_ic50_data.txt",
	sep="\t",
	header=TRUE,
	row.names=1,
	stringsAsFactors=FALSE
	)

# additional information about the cancer cell lines (site etc.)
cell_line_info <- read.table(
	file="gdsc_cell_line_info.txt",
	sep="\t",
	header=TRUE,
	row.names=1,
	stringsAsFactors=FALSE
	)

# find the common set of cell lines in the drug
# and mutation data and sort asciibetically
common_cell_lines <- sort(
	intersect(
		rownames(mutation_data),
		rownames(drug_data)
		)
	)

# make a list holding the tables of data with
# the common cell lines
comb_data <- list(
	mutations = mutation_data[common_cell_lines,],
	drugs = drug_data[common_cell_lines,]
	)

# We can now run and analysis to identify drugs with 
# IC50 values that differ between cell lines with 
# different mutation status or CNV

# try looping over the sets of measurements for
# each mutation/CNV and each drug and then perform
# a spearman rank correlation test
comb_data_all_spearman_results <- data.frame(
	marker=character(),
	drug=character(),
	rho=numeric(),
	pval=numeric(),
	row.names=NULL,
	stringsAsFactors=FALSE
	)
	
r <- 0 # use to add index rows in the data.frame

for(marker in colnames(comb_data$mutations)){
	# marker now has the name of a column containing
	# one of the mutated genes or genes altered by CNA
	for(drug in colnames(comb_data$drugs)){
		# likewise, drug is now the column name of
		# one of the columns in the drugs table
		r <- r + 1
		result <- NULL
		result <- cor.test(
			comb_data$mutations[,marker],
			comb_data$drugs[,drug],
			method="spearman",
			na.action=na.omit
			)
		comb_data_all_spearman_results[r,] <- c(
			marker,
			drug,
			result$estimate,
			result$p.value
			)
	}
}

colnames(comb_data_all_spearman_results) <- c(
	"marker",
	"drug",
	"rho",
	"p.value"
	)

write.table(
	comb_data_all_spearman_results,
	file="comb_data_all_spearman_results.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# now look at the following link for why this is a bad idea
# http://architects.dzone.com/articles/some-tips-r-data-frames

# sort out the types...
comb_data_all_spearman_results <- data.frame(
	marker=as.character(comb_data_all_spearman_results$marker),
	drug=as.character(comb_data_all_spearman_results$drug),
	rho=as.numeric(comb_data_all_spearman_results$rho),
	p.value=as.numeric(comb_data_all_spearman_results$p.value),
	stringsAsFactors=FALSE
	)


# This should be changed for the actual demo but - for now
# the results table has been split into CNV and mutation
# predictors. Read in and visualise


comb_data_mut_spearman_results <- read.table(
	file="comb_data_mut_spearman_results.txt",
	sep="\t",
	header=TRUE,
	stringsAsFactors=FALSE
	)

comb_data_cnv_spearman_results <- read.table(
	file="comb_data_cnv_spearman_results.txt",
	sep="\t",
	header=TRUE,
	stringsAsFactors=FALSE
	)



#
# Visualise the test results
#

# volcano plot
plot(
	comb_data_all_spearman_results$rho,
	-log10(comb_data_all_spearman_results$p.value)
	)


boxplot_drugIC50_by_mutation(
	results=comb_data_mut_spearman_results[which(
		comb_data_mut_spearman_results$p.value <= 0.0000001
		),],
	scores=comb_data$drugs,
	mutations=comb_data$mutations,
	filename="boxplots_drug_IC50_by_mutations.pdf"
	)


boxplot_drugIC50_by_mutation <- function(
	results,
	scores,
	mutations,
	filename
	){
	pdf(file=filename, width=2.5, height=3.5)
	# turn off boxes for plots
	par(bty="n", tcl=-0.2, mai=c(1, 0.95, 0.1, 0.1)) 
	i <- NULL
	for(i in 1:nrow(results)){
		# start by setting all cell lines to wt
		wt_mut_grps_strings <- rep(
			"wt",
			times=length(mutations[,results$marker[i]])
			)
		# set the recurrent/functional mutations
		wt_mut_grps_strings[which(mutations[,results$marker[i]] == 1)] <- "mut"
		wt_grp_rows <- which(wt_mut_grps_strings == "wt")
		func_mut_grp_rows <- which(wt_mut_grps_strings == "mut")
		# boxplot based on all data (wt and mut groups)
		boxplot(
			scores[wt_grp_rows,results$drug[i]],
			scores[func_mut_grp_rows,results$drug[i]],
			pch="",
			names=c("wt", "mutant")
			)
		mtext(results$marker[i], 1, line=2, cex=0.8)
		mtext("status", 1, line=3, cex=0.8)
		mtext(results$drug[i], 2, line=3, cex=0.8)
		mtext("log10 IC50 (µM)", 2, line=2, cex=0.8)
		# jittrered points for each cell line
		if(length(wt_grp_rows) > 0){
			# plot at 1
			points(
				jitter(rep(1,times=length(wt_grp_rows)), amount=0.33),
				scores[wt_grp_rows,results$drug[i]],
				col=rgb(0,0,0,0.25),
				pch=19
				)
		}
		if(length(func_mut_grp_rows) > 0){
			# plot at 2
			points(
				jitter(rep(2,times=length(func_mut_grp_rows)), amount=0.33), scores[func_mut_grp_rows,results$drug[i]], col=rgb(0,0,0,0.25), pch=19
				)
		}
	}		
	dev.off()
}
