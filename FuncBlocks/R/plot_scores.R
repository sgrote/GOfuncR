
# plot the odds-ratio and the 95%-CI from Fishers excact test (two-sided)
# show number of background and candidate genes in pie-charts with size dependent on total nr. of genes

# use as input results from go_enrich and (fwer-threshold or go_ids)

plot_scores = function(res, fwer_threshold=0.05, go_ids=NULL){
	
	root_names = unique(res[[1]][,1])
	root_ids = term[match(root_names, term[,2]) ,4] # TODO: allow custom ontology
	def_root_ids = c("GO:0003674","GO:0005575","GO:0008150") # default root-ids for stable colors
	# TODO:remove default if onto is input

	### check input
	# check that res could be go_enrich-output
	if(!(is.list(res) && all(names(res) == c("results","genes","ref_genome")))){
		stop("Please use an object returned from go_enrich(..., test='wilcoxon') as input (list with 3 elements).")
	}
	# check that it's a wilcox test
	if(colnames(res[[1]][7]) != "FWER_high_rank"){
		stop("Please use the result of a wilcoxon rank sum test performed with go_enrich(..., test='wilcoxon') as input.")
	}
	# check that fwer_threshold is numeric
	if(!is.numeric(fwer_threshold)){
		stop("Please use a numeric fwer_threshold.")
	}

	### define GOs and get annotated genes
	if(is.null(go_ids)){
		# get GOs under FWER-threshold
		go_ids = res[[1]][res[[1]][,"FWER_high_rank"] < fwer_threshold, "node_id"]
		if(length(go_ids) == 0){
			stop(paste("No GO-categories with FWER_high_rank below", fwer_threshold))
		}
	} else {
		# get custom GOs
		go_ids = as.character(go_ids)
	}
	# keep order of input GO's (which gets messed up in get_anno_genes by *apply)
	ordere = data.frame(go_ids, rank=1:length(go_ids))
	
	# get annotated genes for GO-IDs
	in_genes = res[[2]]
	anno = get_anno_genes(go_ids=go_ids, genes=in_genes[,1], ref_genome=res[[3]])
	if(is.null(anno)) return(invisible(anno)) # no annotations - warning from get_anno_genes
	anno$score = in_genes[match(anno$anno_gene, in_genes[,1]), 2]
	# retain original node order
	anno = anno[order(ordere[match(anno$go_id, ordere$go_ids),"rank"]),]

	# get annotated genes for the root nodes
	root_anno = get_anno_genes(go_ids=root_ids, genes=in_genes[,1], ref_genome=res[[3]])
	root_anno$score = in_genes[match(root_anno$anno_gene, in_genes[,1]), 2]

	### Violin Plots
	# colors (some more for possible future custom root-node number option)
	pie_cols = c("#F15A60","#7BC36A","#599BD3","#F9A75B","#9E67AB","#CE7058","#D77FB4")
	root_cols = data.frame(root=def_root_ids, col=pie_cols[1:length(def_root_ids)], stringsAsFactors=FALSE)
	# layout
	layout(matrix(c(1,2),ncol=2), widths=c(5,2))
	ylim = range(c(anno$score, root_anno$score))
	op = par(no.readonly = TRUE) 
	# plot GO-categories
	par(mar=c(6.5,4,4,2), bty="l") #, bty="n") # mar default=c(5, 4, 4, 2)
	violin(anno, root_cols, ylim)
	# plot root nodes
	par(mar=c(6.5,1,4,2), bty="l")
	violin(root_anno, root_cols, ylim, root=TRUE)	
	par(op) 

	out = list(anno_score=anno, anno_score_root=root_anno)
	return(invisible(out))
}


## plot multiple violins side by side
violin = function(plotty, root_cols, ylim, root=FALSE){
	viol_ids = unique(plotty$go_id)
	# find root-node for every GO
	roots = get_names(viol_ids) # also states the root_node name
	roots$root_id = term[match(roots$root_node, term[,2]),4] # add root_id
	# colors
	colore = root_cols[match(roots[,"root_id"],root_cols[,1]),2]
	mainy = ifelse(root, "root nodes", "scores of annotated genes")
	# xlim
	xoff = 0.5
	xlim = c(1-xoff, length(viol_ids)+xoff)
	# plot
	plot(1, type="n", ylim=ylim, xlim=xlim, ylab="score", xaxt="n", xlab="", main=mainy)
	if(root){
		xlabel = roots$root_node
		labelcol = colore
	} else {
		xlabel = viol_ids
		labelcol = "black"
	}
	# x-axis
	xlabpos = ylim[1]-(ylim[2]-ylim[1])/12
	axis(1, at=1:length(viol_ids), labels=FALSE, cex.axis=0.9)
	text(x=1:length(viol_ids), y=xlabpos , labels=xlabel, srt=45, adj=c(1,1), xpd=TRUE, cex=0.8, col=labelcol)
	for(i in 1:length(viol_ids)){
		scores = plotty[plotty$go_id == viol_ids[i], "score"]
		if (length(scores)==1) {
			points(x=i, y=scores, col=colore[i], pch=16, cex=1.5)
		} else {
			vioplot(scores, at=i, col=colore[i], add=TRUE)
		}
	}
	# n genes
	n = table(plotty$go_id)
	text(x=1:length(viol_ids), y=ylim[2], labels=paste0("n=",n[viol_ids]), xpd=TRUE, cex=0.8, col=labelcol, pos=3, offset=0.5)	
}
