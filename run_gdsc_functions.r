# ============================================ #
# ICR BIC GDSC workshop
# Functions...
# https://github.com/DrJCampbell/genexdrugdemo
# ============================================ #


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
				jitter(rep(2,times=length(func_mut_grp_rows)), amount=0.33),
				scores[func_mut_grp_rows,results$drug[i]],
				col=rgb(0,0,0,0.25),
				pch=19
				)
		}
	}		
	dev.off()
}


boxplot_drugIC50_by_copy_number <- function(
	results,
	scores,
	mutations,
	filename
	){
	pdf(file=filename, width=10, height=3.5)
	# turn off boxes for plots
	par(bty="n", tcl=-0.2, mai=c(1, 0.95, 0.1, 0.1)) 
	i <- NULL
	for(i in 1:nrow(results)){
		# boxplot IC50 dependent on CP
		boxplot(
			scores[,results$drug[i]] ~ mutations[,results$marker[i]]
			)
		mtext(results$marker[i], 1, line=2, cex=0.8)
		mtext("copy number", 1, line=3, cex=0.8)
		mtext(results$drug[i], 2, line=3, cex=0.8)
		mtext("log10 IC50 (µM)", 2, line=2, cex=0.8)
	}		
	dev.off()
}


scatterplot_drugIC50_by_copy_number <- function(
	results,
	scores,
	mutations,
	filename
	){
	pdf(file=filename, width=10, height=3.5)
	# turn off boxes for plots
	par(bty="n", tcl=-0.2, mai=c(1, 0.95, 0.1, 0.1)) 
	i <- NULL
	for(i in 1:nrow(results)){
		# boxplot IC50 dependent on CP
		plot(
			scores[,results$drug[i]] ~ mutations[,results$marker[i]],
			pch=19,
			col=rgb(0,0,0,0.25),
			xlab=paste(results$marker[i], "copy number"),
			ylab=paste(results$drug[i], "IC50")
			)
#		mtext(results$marker[i], 1, line=2, cex=0.8)
#		mtext("copy number", 1, line=3, cex=0.8)
#		mtext(results$drug[i], 2, line=3, cex=0.8)
#		mtext("log10 IC50 (µM)", 2, line=2, cex=0.8)
	}		
	dev.off()
}