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


