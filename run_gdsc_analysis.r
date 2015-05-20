#
# ICR BIC GDSC workshop
#

# can change the next line to the
# directory where your files are
setwd("/Users/jamesc/Dropbox/Bioinformatics_course_materials")




# read in the expression, mutation
# and copy number data for cell lines
#predictors <- read.table(
#	file="en_input_w5.csv",
#	header=TRUE,
#	sep=",",
#	row.names=1
#	)


mutation_data <- read.table(
	file="gdsc_CNV_mutation_predictors.txt",
	sep="\t",
	header=TRUE,
	row.names=1
	)

drug_data <- read.table(
	file="gdsc_ic50_data.txt",
	sep="\t",
	header=TRUE,
	row.names=1,
	stringsAsFactors=FALSE
	)

cell_line_info <- read.table(
	file="gdsc_cell_line_info.txt",
	sep="\t",
	header=TRUE,
	row.names=1,
	stringsAsFactors=FALSE
	)


