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
	stringsAsFactors=FALSE
	)
colnames(comb_data_all_spearman_results) <- c(
	"marker",
	"drug",
	"rho",
	"p.value"
	)
for(marker in colnames(comb_data$mutations)){
	# marker now has the name of a column containing
	# one of the mutated genes or genes altered by CNA
	for(drug in colnames(comb_data$drugs)){
		# likewise, drug is now the column name of
		# one of the columns in the drugs table
		result <- NULL
		result <- cor.test(
			comb_data$mutations[,marker],
			comb_data$drugs[,drug],
			method="spearman",
			na.action=na.omit
			)
		comb_data_all_spearman_results <- rbind(
			comb_data_all_spearman_results,
			cbind(
				marker,
				drug,
				result$estimate,
				result$p.value
				)
			)
	}
}













