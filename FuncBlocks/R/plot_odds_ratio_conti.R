
# plot the odds-ratio and the 95%-CI from Fishers excact test (two-sided)
# show number of background and candidate genes in pie-charts with size dependent on total nr. of genes

# use as input results from go_enrich and (fwer-threshold or go_ids)

# TODO: test plotted odds ratio with an example

plot_odds_ratio_conti = function(res, fwer_threshold=0.05, go_ids=NULL){
	
	### check input
	# check that res could be go_enrich-output
	if(!(is.list(res) && all(names(res) == c("results","genes","ref_genome")))){
		stop("Please use an object returned from go_enrich as input (list with 3 elements).")
	}
	
	
	# check that fwer_threshold is numeric
	if(!is.numeric(fwer_threshold)){
		stop("Please use a numeric fwer_threshold.")
	}
	
	### define GOs and get annotated genes
	if(is.null(go_ids)){
		# get GOs under FWER-threshold
		go_ids = res[[1]][res[[1]][,7] < fwer_threshold, "node_id"]
		if(length(go_ids) == 0){
			stop(paste("No GO-categories with FWER_overrep below", fwer_threshold))
		}
	} else {
		# get custom GOs
		go_ids = as.character(go_ids)
	}
	# restrict to input genes
	in_genes = res[[2]]
	anno = get_anno_genes(go_ids=go_ids, genes=in_genes[,1], ref_genome=res[[3]])
	if(is.null(anno)) return(invisible(anno)) # no annotations - warning from get_anno_genes
	# add scores
	anno = cbind(anno, in_genes[match(anno[,2], in_genes[,1]), 2:ncol(in_genes)])
	# sum scores of all genes for one GO
	anno = aggregate(anno[,3:6], list(go_id=anno[,1]), sum)

	# keep order of input GO's (which gets messed up in get_anno_genes by *apply)
	ordere = data.frame(go_ids, rank=1:length(go_ids))

	# find root-node for every GO
#	go_ids = unique(anno[,1]) # remove any GO without annotations # TODO: remove?
	
	### perform fishers exact test for every GO
	fish_odds = data.frame(anno, t(apply(anno, 1, fisher_conti)))
	colnames(fish_odds)[c(2:9)] = c("A","B","C","D","odds_ratio","ci95_low","ci95_high","p")

	# recreate original order of Go-ids # TODO: needed?
	fish_odds = fish_odds[order(ordere[match(fish_odds$go_id, ordere$go_ids),"rank"]),]
	out = fish_odds

	# replace Inf or -Inf CI for plotting
	inf_go = fish_odds[is.infinite(fish_odds$odds_ratio),"go_id"]
	if(length(inf_go) > 0){
		warning(paste("The following GOs have an infinite odds ratio due to no background gene annotation: "), paste(inf_go, collapse=", "), sep="")
	}
	zero_go = fish_odds[fish_odds$odds_ratio==0,"go_id"]
	if(length(zero_go) > 0){
		warning(paste("The following GOs have an odds ratio of 0: "), paste(zero_go, collapse=", "), sep="")
	}
	############# TODO influenced by setting par(oma) before function call?
	
	### plot the odds-ratios, with CI
	op = par(no.readonly = TRUE) 
	# 3 panels
	layout(matrix(c(1:2),ncol=1), widths=c(5),heights=c(3,3))

	# odds-ratio
	par(mar=c(0.5,4,3,2), bty="l") #, bty="n")
	ymin = min(c(1,fish_odds$ci95_low[is.finite(fish_odds$ci95_low) & fish_odds$odds_ratio!=0]))
	ymax = max(c(1,fish_odds$ci95_high[is.finite(fish_odds$ci95_high) & fish_odds$odds_ratio!=0]))
	suppressWarnings(plot(fish_odds$odds_ratio, pch=19, ylab="", xaxt="n", xlab="", main="odds ratio (A/B) / (C/D)", xlim=c(0.5,nrow(fish_odds)+1), ylim=c(ymin,ymax), cex.axis=0.8, log="y", las=2, 
	panel.first={grid(0, NULL, lty=1, col=colors()[2])}))
	# 95%-CI
	suppressWarnings(arrows(c(1:nrow(fish_odds),1:nrow(fish_odds)),c(fish_odds$ci95_high,fish_odds$ci95_low), c(1:nrow(fish_odds),1:nrow(fish_odds)), c(fish_odds$odds_ratio, fish_odds$odds_ratio), angle=90, code=1, length=0.03))
	# horizontal line at 1
	abline(h=1, col="#F15A60")
	axis(4,at=1,labels="1", col="#F15A60", col.axis="#F15A60", las=1, cex.axis=0.8)
	
	# pie charts
	pie_cols = colors()[c(124,132,59,137)]
#	pie_cols = colors()[c(592,566,503,556)]

	par(mar=c(5.5,4,3,2), bty="l")
	plot(1,xlim=c(0.5,nrow(fish_odds)+1), ylim=c(0,2),type="n", main="annotated genes", xlab="", ylab="", xaxt="n", yaxt="n")
	radi_units = 0.4/log(max(rowSums(fish_odds[,2:5]))+1)
	for(i in 1:nrow(fish_odds)){
		add.pie(z=as.numeric(fish_odds[i,2:3]), x=i, y=1.6, radius=log(sum(fish_odds[i,2:3])+1)*radi_units, labels="", col=pie_cols[1:2])
		add.pie(z=as.numeric(fish_odds[i,4:5]), x=i, y=0.6, radius=log(sum(fish_odds[i,4:5])+1)*radi_units, labels="", col=pie_cols[3:4])
		
	}
	legend("right", fill=pie_cols, bty="n", legend=c("A","B","C","D"))
	text(x=1:nrow(fish_odds), y=1.9, labels=paste0("n=",rowSums(fish_odds[,2:5])), xpd=TRUE, cex=0.6, pos=3, offset=0.5)
	axis(1, at=1:nrow(fish_odds), labels=FALSE, cex.axis=0.8)
	text(x=1:nrow(fish_odds), y=-0.35, labels=fish_odds$go_id, srt=45, adj=1, xpd=TRUE, cex=0.8)
	
	par(op)
	return(invisible(out))
}


fisher_conti = function(fish_line){
	# contingency table
	conti = matrix(as.numeric(fish_line[2:5]), ncol=2)
	fish = fisher.test(conti)
	out = c(fish$estimate, fish$conf.int, fish$p.value) 
	return(out)	
}

