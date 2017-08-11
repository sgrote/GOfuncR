
# plot the p(A/B) and the 95%-CI from Fishers excact test (two-sided)
# show number of background and candidate genes in pie-charts with size dependent on total nr. of genes

# use as input results from go_enrich and (fwer-threshold or go_ids)

plot_binom = function(res, fwer_threshold=0.05, go_ids=NULL){
	

	### check input
	# check that res could be go_enrich-output
	if(!(is.list(res) && all(names(res) == c("results","genes","ref_genome")))){
		stop("Please use an object returned from go_enrich as input (list with 3 elements).")
	}
	# check that fwer_threshold is numeric
	if(!is.numeric(fwer_threshold)){
		stop("Please use a numeric fwer_threshold.")
	}
	# get IDs for root_nodes from res
	root_names = unique(res[[1]][,1])
	root_ids = term[match(root_names, term[,2]) ,4] # TODO: allow custom ontology
	def_root_ids = c("GO:0003674","GO:0005575","GO:0008150") # default root-ids for stable colors
	# TODO:remove default if onto is input
	
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

	# background defined: restrict to input genes
	# restrict to input genes
	in_genes = res[[2]]
	anno = get_anno_genes(go_ids=go_ids, genes=in_genes[,1], ref_genome=res[[3]])
	if(is.null(anno)) return(invisible(anno)) # no annotations - warning from get_anno_genes
	# add scores
	anno = cbind(anno, in_genes[match(anno[,2], in_genes[,1]), 2:ncol(in_genes)])
	# sum scores of all genes for one GO
	anno = aggregate(anno[,3:4], list(go_id=anno[,1]), sum) # TODO: ncol
	
	# keep order of input GO's (which gets messed up in get_anno_genes by *apply)
	ordere = data.frame(go_ids, rank=1:length(go_ids))

	# find root-node for every GO
	go_ids = unique(anno[,1]) # remove any GO without annotations
	roots = get_names(go_ids) # also states the root_node
	roots$root_id = term[match(roots$root_node, term[,2]),4]

	# get annotated genes for the root nodes
	#TODO: allow different root nodes (term has "all" root...,remove that) -> 
	# -> require root nodes-vector as input from user with default values; or use only roots of required GOs? 
	root_anno = get_anno_genes(go_ids=root_ids, genes=in_genes[,1], ref_genome=res[[3]])
	root_anno = cbind(root_anno, in_genes[match(root_anno[,2], in_genes[,1]), 2:ncol(in_genes)])
	# compute p_A_root for root-nodes
	root_anno = aggregate(root_anno[,3:4], list(go_id=root_anno[,1]), sum) # TODO: ncol
	root_anno$p_A_root = root_anno[,2] / (root_anno[,2]+root_anno[,3])
		
	# perform binomial test for every GO
	anno$root_id = roots[match(anno$go_id, roots$go_id), "root_id"]
	anno$p_A_root = root_anno[match(anno$root_id, root_anno$go_id), "p_A_root"]
	binom = data.frame(anno, t(apply(anno[,c(2,3,5)], 1, binomy)))
	colnames(binom)[c(6:9)] = c("p_A_node","ci95_low","ci95_high","p")

	# recreate original order of Go-ids
	binom = binom[order(ordere[match(binom$go_id, ordere$go_ids),"rank"]),]
	out = binom

	############# TODO influenced by setting par(oma) before function call?
	
	### plot p_A_node, with CI
	op = par(no.readonly = TRUE) 
	# 3 panels
	layout(matrix(c(1,2,3,3),ncol=2), widths=c(4,1),heights=c(3,2))

	# p_A_node
	par(mar=c(0.5,4,3,1), bty="l") #, bty="n")
	rangy = range(c(binom$ci95_low, binom$ci95_high, root_anno$p_A_root))
	ymin = max(0, rangy[1] - 0.2 * (rangy[2] - rangy[1]))
	ymax = min(1, rangy[2] + 0.2 * (rangy[2] - rangy[1]))
#	ymin = 0
#	ymax = 1
	suppressWarnings(plot(binom$p_A_node, pch=19, ylab="", xaxt="n", xlab="", main="p(A) in node", xlim=c(0.5,nrow(binom)+0.5), ylim=c(ymin,ymax), cex.axis=0.8, las=2, 
	panel.first={grid(0, NULL, lty=1, col=colors()[2])}))
	
	# 95%-CI
	suppressWarnings(arrows(c(1:nrow(binom),1:nrow(binom)),c(binom$ci95_high,binom$ci95_low), c(1:nrow(binom),1:nrow(binom)), c(binom$p_A_node, binom$p_A_node), angle=90, code=1, length=0.03))

	# horizontal line for root-node p_A
	pie_cols = c("#F15A60","#7BC36A","#599BD3","#F9A75B","#9E67AB","#CE7058","#D77FB4")
	root_cols = data.frame(root=def_root_ids,col=pie_cols[1:length(def_root_ids)],stringsAsFactors=FALSE)
	abline(h=root_anno$p_A_root, col=root_cols[match(root_anno$go_id, root_cols$root),"col"])

	
	# pie charts
	par(mar=c(5.5,4,3,1), bty="l")

	plot(1,xlim=c(0.5,nrow(binom)+0.5), ylim=c(0,1),type="n", main="annotated scores", xlab="", ylab="", xaxt="n", yaxt="n")
	radi_units = 0.4/log(max(rowSums(binom[,2:3]))+1)
	pie_cols = colors()[c(124,132,59,137)]
	for(i in 1:nrow(binom)){
		a = binom[i,2]
		b = binom[i,3]
		add.pie(z=c(b,a), x=i, y=0.6, radius=log(a+b+1)*radi_units, labels="", col=c("#737373",root_cols[match(binom[i,"root_id"],root_cols[,1]),2]))
	}
	text(x=1:nrow(binom), y=0.08, labels=paste(binom[,2],rowSums(binom[,2:3]),sep=" / "),col= root_cols[match(binom[,"root_id"],root_cols[,1]),2], xpd=TRUE, cex=0.8)
	axis(1, at=1:nrow(binom), labels=FALSE, cex.axis=0.8)
	text(x=1:nrow(binom), y=-0.25, labels=binom$go_id, srt=45, adj=1, xpd=TRUE, cex=0.8)
	
	# pie charts for root nodes
	par(mar=c(5.5,1,3,1), bty="l")
	plot(1,xlim=c(-0.3,2), ylim=c(0.5,3.5),type="n", main="annotated genes\nroot nodes", xlab="", ylab="", xaxt="n", yaxt="n", yaxs="i")
	radi_units = 0.3/log(max(rowSums(root_anno[,2:3]))+1)
	for(i in 1:nrow(root_anno)){
		a = root_anno[i,2]
		b = root_anno[i,3]
		add.pie(z=c(b,a), x=1, y=i, radius=log(b+a+1)*radi_units, labels="", col=c("#737373",root_cols[match(root_anno[i,1],root_cols[,1]),2]))
	}
	text(x=1, y=(0.4 + 1:nrow(root_anno)), labels=get_names(root_anno$go_id)[,2], col=root_cols[match(root_anno[,1],root_cols[,1]),2], cex=0.8)
	text(x=1, y=(-0.4 + 1:nrow(root_anno)), labels=paste(root_anno[,2],rowSums(root_anno[,2:3]), sep=" / "), cex=0.8, col=root_cols[match(root_anno[,1],root_cols[,1]),2])
	
	par(op)
	# TODO: also return p? add significance asterix for fisher-test? - but does not account for multiple testing, # rather FWER-over and under?
	return(invisible(out))
}


binomy = function(bino_line){
	# binom.test(c(success, failure), p)
	bino_line = unname(bino_line)
	bino = binom.test(c(bino_line[1], bino_line[2]), p=bino_line[3])
	out = c(bino$estimate, bino$conf.int, bino$p.value) 
	return(out)	
}

