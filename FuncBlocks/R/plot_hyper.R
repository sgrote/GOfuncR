
# plot the odds-ratio and the 95%-CI from Fishers excact test (two-sided)
# show number of background and candidate genes in pie-charts with size dependent on total nr. of genes
# take table with candidate and background gene numbers prepared in plot_anno_scores

plot_hyper = function(aggrego, root_aggrego){
	
	
	print(head(aggrego))
	print(head(root_aggrego))
	
	# flip canidate and background columns
	aggrego[,c(2:3, 5:6)] = aggrego[,c(3,2,6,5)]
	colnames(aggrego)[c(2:3, 5:6)] = c("candi_genes", "bg_genes", "root_candi_genes", "root_bg_genes")
	root_aggrego[,c(2:3)] = root_aggrego[,c(3,2)]
	colnames(root_aggrego)[c(2:3)] = c("candi_genes", "bg_genes")
	
	print(head(aggrego))
	print(head(root_aggrego))
	
	fish_odds = data.frame(aggrego, t(apply(aggrego, 1, fisher)))
	print(fish_odds)
	out = fish_odds
#	fish_odds = fish_odds[,c(1,2,6,5,4,3,8,7,10:13)]
#	colnames(fish_odds)[c(2,6,9:12)] = c("go_name","root_name","odds_ratio","ci95_low","ci95_high","p")
#	# recreate original order of Go-ids
#	fish_odds = fish_odds[order(ordere[match(fish_odds$go_id, ordere$go_ids),"rank"]),]
#	out = fish_odds

#	# replace Inf or -Inf CI for plotting
#	inf_go = fish_odds[is.infinite(fish_odds$odds_ratio),"go_id"]
#	if(length(inf_go) > 0){
#		warning(paste("The following GO-categories have an infinite odds ratio due to no background gene annotation: "), paste(inf_go, collapse=", "), sep="")
#	}
#	zero_go = fish_odds[fish_odds$odds_ratio==0,"go_id"]
#	if(length(zero_go) > 0){
#		warning(paste("The following GOs have an odds ratio of 0: "), paste(zero_go, collapse=", "), sep="")
#	}
#	############# TODO influenced by setting par(oma) before function call?
	
#	### plot the odds-ratios, with CI
#	op = par(no.readonly = TRUE) 
#	# 3 panels
#	layout(matrix(c(1,2,3,3),ncol=2), widths=c(4,1),heights=c(3,2))

#	# odds-ratio
#	par(mar=c(0.5,4,3,1), bty="l") #, bty="n")
#	ymin = min(c(1,fish_odds$ci95_low[is.finite(fish_odds$ci95_low) & fish_odds$odds_ratio!=0]))
#	ymax = max(c(1,fish_odds$ci95_high[is.finite(fish_odds$ci95_high) & fish_odds$odds_ratio!=0]))
#	suppressWarnings(plot(fish_odds$odds_ratio, pch=19, ylab="", xaxt="n", xlab="", main="odds ratio", xlim=c(0.5,nrow(fish_odds)+0.5), ylim=c(ymin,ymax), cex.axis=0.8, log="y", las=2, 
#	panel.first={grid(0, NULL, lty=1, col=colors()[2])}))
#	# 95%-CI
#	suppressWarnings(arrows(c(1:nrow(fish_odds),1:nrow(fish_odds)),c(fish_odds$ci95_high,fish_odds$ci95_low), c(1:nrow(fish_odds),1:nrow(fish_odds)), c(fish_odds$odds_ratio, fish_odds$odds_ratio), angle=90, code=1, length=0.03))
#	# horizontal line at 1
#	abline(h=1, col="#F15A60")
#	axis(4,at=1,labels="1", col="#F15A60", col.axis="#F15A60", las=1, cex.axis=0.8)
	
#	# pie charts
#	par(mar=c(5.5,4,3,1), bty="l")
#	pie_cols = c("#F15A60","#7BC36A","#599BD3","#F9A75B","#9E67AB","#CE7058","#D77FB4")
#	root_cols = data.frame(root=def_root_ids,col=pie_cols[1:length(def_root_ids)],stringsAsFactors=FALSE)

#	plot(1,xlim=c(0.5,nrow(fish_odds)+0.5), ylim=c(0,1),type="n", main="annotated genes", xlab="", ylab="", xaxt="n", yaxt="n")
#	radi_units = 0.4/log(max(rowSums(fish_odds[,3:4]))+1)
#	for(i in 1:nrow(fish_odds)){
#		w = fish_odds[i,3]
#		b = fish_odds[i,4]
#		add.pie(z=c(b,w), x=i, y=0.6, radius=log(b+w+1)*radi_units, labels="", col=c("#737373",root_cols[match(fish_odds[i,"root_id"],root_cols[,1]),2]))
#	}
#	text(x=1:nrow(fish_odds), y=0.08, labels=paste(fish_odds[,3],rowSums(fish_odds[,3:4]),sep=" / "),col= root_cols[match(fish_odds[,"root_id"],root_cols[,1]),2], xpd=TRUE, cex=0.8)
#	axis(1, at=1:nrow(fish_odds), labels=FALSE, cex.axis=0.8)
#	text(x=1:nrow(fish_odds), y=-0.25, labels=fish_odds$go_id, srt=45, adj=1, xpd=TRUE, cex=0.8)
	
#	# pie charts for root nodes
#	par(mar=c(5.5,1,3,1), bty="l") # TODO, replace 3 with number of root nodes
#	plot(1,xlim=c(-0.3,2), ylim=c(0.5,3.5),type="n", main="annotated genes\nroot nodes", xlab="", ylab="", xaxt="n", yaxt="n", yaxs="i")
#	radi_units = 0.3/log(max(rowSums(urn_table[,1:2]))+1)
#	for(i in 1:nrow(urn_table)){
#		b = urn_table[i,1]
#		w = urn_table[i,2]
#		add.pie(z=c(b,w), x=1, y=i, radius=log(b+w+1)*radi_units, labels="", col=c("#737373",root_cols[match(rownames(urn_table)[i],root_cols[,1]),2]))
#	}
#	text(x=1, y=(0.4 + 1:nrow(urn_table)), labels=urn_table$name, col=root_cols[match(rownames(urn_table),root_cols[,1]),2], cex=0.8)
#	text(x=1, y=(-0.4 + 1:nrow(urn_table)), labels=paste(urn_table[,2],rowSums(urn_table[,1:2]), sep=" / "), cex=0.8, col=root_cols[match(rownames(urn_table),root_cols[,1]),2])
	
#	par(op)
	
#	 #TODO: also return p? add significance asterix for fisher-test? - but does not account for multiple testing, # rather FWER-over and under?
#	return(invisible(out))
}


fisher = function(fish_line){
	# drawn white and black balls in node
	w = as.numeric(fish_line["candi_genes"])
	b = as.numeric(fish_line["bg_genes"])
	# white and black not drawn (left in urn)
	wn = as.numeric(fish_line["root_candi_genes"]) - w
	bn = as.numeric(fish_line["root_bg_genes"]) - b
	# contingency table
	conti = cbind(c(w,b),c(wn,bn))
	fish = fisher.test(conti)
	out = c(fish$estimate, fish$conf.int, fish$p.value) 
	return(out)	
}

