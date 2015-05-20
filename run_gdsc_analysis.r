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


# volcano plot
plot(
	comb_data_all_spearman_results$rho,
	-log10(comb_data_all_spearman_results$p.value)
	)

# which are the strongest associations (drug sensitivities)
# r ≤ -0.1, p ≤ 0.05...
comb_data_all_spearman_results[
	which(
		comb_data_all_spearman_results$rho >= 0.1 &
		comb_data_all_spearman_results$p.value <= 0.001
		)
		,
	]

write.table(
	comb_data_all_spearman_results,
	file="comb_data_all_spearman_results.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)









