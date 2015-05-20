genexdrugdemo
=============
Demonstration R code to test for differential sensitivity of cancer cell lines to drugs, explained by genetic information


Running the demo
----------------
The R code in the file named 'run_gdsc_analysis.r'. Open R-Studio or your favourite R environment and run the code in blocks to see what happens. You will need to modify the setwd() line so that you read the files from the correct place (the folder where this file is).

Getting the data
----------------
The data were downloaded from the ftp.sanger.ac.uk FTP server using curl

	curl \
	ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_en_input_w5.csv.zip \
	> gdsc_en_input_w5.csv.zip

Some clean-up / re-formatting of the files was done to make life easier in R. Roughly, the following steps were taken

+ strip all non-alphanumeric characters from targets header row
+ strip _IC_50 suffix from all compounds in header row
+ concatenate compounds and targets to make the actual names used for R analysis
+ beware line endings - everything was converted to Unix (\n) from Mac (\r)
+ the mutation data was filtered to remove all expression data
+ the remaining mutation and CNV data was transposed so that variables became columns
+ probably lots of other things that were not well documented
