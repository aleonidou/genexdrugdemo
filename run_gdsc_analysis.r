# ============================================ #
# ICR BIC GDSC workshop
# https://github.com/DrJCampbell/genexdrugdemo
# ============================================ #

# change the next line to the directory where your files are
setwd("./")

# run the next line if gplots is not installed
# install.packages("gplots", dependencies=TRUE)
require(gplots)

# source (read in) some code from the following file
source("./run_gdsc_functions.r")

# read a file with gene mutation and copy number data
mutation_data <- read.table(
	file="gdsc_CNV_mutation_predictors.txt",
	sep="\t",
	header=TRUE,
	row.names=1
	)
	# mutation data now contains a table (data.frame) with
	# 624 rows corresponding to cancer cell lines and 510
	# columns corresponding to either mutation status [0/1]
	# or infered absolute copy number values

# read file with drug IC50 values
drug_data <- read.table(
	file="gdsc_ic50_data.txt",
	sep="\t",
	header=TRUE,
	row.names=1,
	stringsAsFactors=FALSE
	)
	# Similarly, mutation data now contains a table 
	# 707 rows corresponding to cancer cell lines and 140
	# columns corresponding to IC50 values of compounds



# find the common set of cell lines in the drug
# and mutation data and sort asciibetically
common_cell_lines <- sort(
	intersect(
		rownames(mutation_data),
		rownames(drug_data)
		)
	)
	# sort() sorts a list of things
	# intersect() finds the common things in two lists
	# rownames() gets the names of the rows in a table

# make a list holding the tables of data with
# the common cell lines
comb_data <- list(
	mutations = mutation_data[common_cell_lines,],
	drugs = drug_data[common_cell_lines,]
	)
	# mutation_data[common_cell_lines,] chooses rows
	# from mutation_data that are in common_cell_lines
	# and orders the rows as per common_cell_lines


#
# visualise the data sets
#

# create a colour pallete and breaks (data bin ranges)
ic50_breaks=seq(-10, 10, by=0.2) 
ic50_breaks=append(ic50_breaks, 40)
ic50_breaks=append(ic50_breaks, -40, 0)
ic50_col <- colorpanel(n=length(ic50_breaks)-1,low="blue", high="white")

# Make a really basic heatmap of the drug IC50s
pdf("ic50_data_heatmap_default.pdf", width=50, height=50)
image(
	as.matrix(comb_data$drugs),
	col=ic50_col
	)
dev.off()

# Make a heatmap of mutation and copy number values
cnv_breaks=seq(0, 10, by=1) 
cnv_col <- colorpanel(n=length(cnv_breaks)-1,low="white",  high="red")
pdf("mutation_data_heatmap_default.pdf", width=50, height=50)
heatmap(
	as.matrix(comb_data$mutations),
	breaks=cnv_breaks,
	col=cnv_col
	)
dev.off()



# We can now run an analysis to identify drugs with 
# IC50 values that differ between cell lines with 
# different mutation status or CNV

# look at the distribution
hist(
	comb_data$drugs[,"PLX4720_RAF"],
	breaks=50,
	xlab="PLX4720 (RAFi) IC50 (µM)",
	main="Histogram showing PLX4720 log10 IC50 distribution"
	)

# non-parametric correlation test
cor_example <- cor.test(
	comb_data$mutations[,"BRAF.MUT"],
	comb_data$drugs[,"PLX4720_RAF"],
	na.rm=TRUE,
	method="spearman"
	)

# look at what is available in cor_example
str(cor_example)

# output visualisations of RAFi response dependent on BRAF mutation to a PDF
pdf("box_and_stripchart_BRAF_RAFi.pdf", width=3, height=4.5)
	# visualise the distribution by group using stripchart
	stripchart(
		comb_data$drugs[,"PLX4720_RAF"] ~ comb_data$mutations[,"BRAF.MUT"],
		vertical=TRUE,
		pch=19,
		col=rgb(0,0,0,0.5),
		method="jitter",
		group.names=c("WT", "mutant"),
		ylab="PLX4720 (RAFi) log10 IC50 (µM)"
		)

	# visualise the distribution by group using boxplot
	boxplot(
		comb_data$drugs[,"PLX4720_RAF"] ~ comb_data$mutations[,"BRAF.MUT"],
		group.names=c("WT", "mutant"),
		ylab="PLX4720 (RAFi) log10 IC50 (µM)"
		)
dev.off()

# ggplot module is good for highly customisable data visualisation


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
		
		r <- r + 1		# initially r goes from
							# 0 to 1 (first row), then
							# r increments by 1
		
		result <- NULL	# clear any previous result
								# in case cor.test() fails
								# and the result is stale
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

# now look at the following link for why this is not considered
# good R coding style. 
# http://architects.dzone.com/articles/some-tips-r-data-frames

# sort out the types... marker and drug are character strings
# rho and p.value are numeric
comb_data_all_spearman_results <- data.frame(
	marker=as.character(comb_data_all_spearman_results$marker),
	drug=as.character(comb_data_all_spearman_results$drug),
	rho=as.numeric(comb_data_all_spearman_results$rho),
	p.value=as.numeric(comb_data_all_spearman_results$p.value),
	stringsAsFactors=FALSE
	)

# or just read the file back in as a dataframe
comb_data_all_spearman_results <- read.table(
	file="comb_data_all_spearman_results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

# separate the predictor (marker) rows the end in
# '.CP' from those that end in '.MUT' and any others

comb_data_cnv_spearman_results <- comb_data_all_spearman_results[
	grep(".CP", comb_data_all_spearman_results$marker),]

comb_data_mut_spearman_results <- comb_data_all_spearman_results[
	grep(".MUT", comb_data_all_spearman_results$marker),]


# To save time for the demo, the correlation tests have
# already been run. We can read these in and visualise the 
# results


# read in the previously run test results for the
# mutation data
comb_data_mut_spearman_results <- read.table(
	file="comb_data_mut_spearman_results.txt",
	sep="\t",
	header=TRUE,
	stringsAsFactors=FALSE
	)

# read in the previously run test results for the
# copy number data
comb_data_cnv_spearman_results <- read.table(
	file="comb_data_cnv_spearman_results.txt",
	sep="\t",
	header=TRUE,
	stringsAsFactors=FALSE
	)

#
# Visualise the test results
#


# call a function (defined in the file we sourced earlier)
# that ouputs a pdf file with multiple boxplots showing
# IC50 for drugs in mutant and non-mutant groups of cell
# lines
boxplot_drugIC50_by_mutation(
	results=comb_data_mut_spearman_results[which(
		comb_data_mut_spearman_results$p.value <= 1e-06
		),],
	scores=comb_data$drugs,
	mutations=comb_data$mutations,
	filename="boxplots_drug_IC50_by_mutations.pdf"
	)

# call another function that makes boxplots of
# drug IC50 where cell lines were split by gene
# copy number
boxplot_drugIC50_by_copy_number(
	results=comb_data_cnv_spearman_results[which(
		comb_data_cnv_spearman_results$p.value <= 0.001
		),],
	scores=comb_data$drugs,
	mutations=comb_data$mutations,
	filename="boxplots_drug_IC50_by_copy_number.pdf"
	)


# call another function that makes scatter plots
# of drug IC50 where cell lines were split by gene
# copy number
scatterplot_drugIC50_by_copy_number(
	results=comb_data_cnv_spearman_results[which(
		comb_data_cnv_spearman_results$p.value <= 0.001
		),
	],
	scores=comb_data$drugs,
	mutations=comb_data$mutations,
	filename="scatterplots_drug_IC50_by_copy_number.pdf"
	)


# ======= #
# The end
# ======= #
