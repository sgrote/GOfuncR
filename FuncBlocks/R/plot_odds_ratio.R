
# plot the odds-ratio and the 95%-CI from Fishers excact test (two-sided)
# show number of background and candidate genes in pie-charts with size dependet on total nr. of genes

# use as input results from go_enrich and (fwer-threshold or go_ids)

plot_odds_ratio = function(res, fwer_threshold=0.05, go_ids=NULL){

	### check input
	# check that res could be go_enrich-output
	if(!(is.list(res) && all(names(res) == c("results","genes")))){
		stop("Please use an object returned from go_enrich as input (list with 2 elements).\n Alternatively go_ids need to be defined.")
	}
	# check that it's a hypergeometric test
	in_genes = res[[2]]
	if(!(all(in_genes %in% c(1,0)))){
		stop("Please use the result of an hypergeometric test performed with go_enrich as input.")
	}
	# check if background is defined
	bgdef = TRUE
	if(all(in_genes==1)){
		bgdef = FALSE
	}
	
	### define GOs and get annotated genes
	if(is.null(go_ids)){
		# get GOs under FWER-threshold
		anno = get_anno_genes(res, fwer_threshold, background=TRUE)
		anno = anno[,-3]
	} else {
		# get custom GOs, add score 1/0
		go_ids = as.character(go_ids)
		if (bgdef){
			# background defined: restrict to input genes
			anno = get_anno_genes(go_ids=go_ids, genes=names(in_genes))
		} else {
			# background not defined: get all genes annotations
			anno = get_anno_genes(go_ids=go_ids)
		}	
		anno$score = 0
		anno[anno[,2] %in% names(in_genes[in_genes==1]), "score"] = 1
	}
	
	# find root-node for every GO
	go_ids = unique(anno[,1]) # if is.null(go_ids) define GOs; else remove any GO without annotations
	roots = get_names(go_ids) # also states the root_node
	roots$root_id = term[match(roots$root_node, term[,2]),4]

	# get annotated genes for the root nodes
	root_anno = get_anno_genes(go_ids=unique(roots$root_id))
	root_anno$score = 0
	root_anno[root_anno[,2] %in% names(in_genes[in_genes==1]), "score"] = 1
	
	### perform fishers exact test for every GO
	drawn_table = do.call(rbind, tapply(anno$score, anno$go_id, table))
	colnames(drawn_table) = c("b","w")
	drawn_table = data.frame(roots, drawn_table[match(roots$go_id, rownames(drawn_table),),])
	urn_table = do.call(rbind, tapply(root_anno$score, root_anno$go_id, table))
	fish_odds = data.frame(drawn_table[,1], t(apply(drawn_table, 1, fisher)))
	colnames(fish_odds) = c("go_id", "odds_ratio", "CI_low", "CI_high","p")
	
	### plot the odds-ratios, with CI
	ymin = min(0,log(min(fish_odds[,3]),10))
	ymax = log(max(fish_odds[,4]),10)
	# y-range for upper pane
	rup = ymax-ymin
	# y-range for entire plot	
	rangy = c(ymin-(rup/2), ymax)
	# odds-ratio
	plot(log(fish_odds$odds_ratio,10), ylim=rangy, pch=19, yaxt="n", ylab="", xaxt="n", xlab="", main="odds-ratio", xlim=c(0.5,nrow(fish_odds)+0.5), cex.axis=0.8)
	# x-axis
	axis(1, at=1:nrow(fish_odds), labels=FALSE, cex.axis=0.8)
	text(x=1:nrow(fish_odds), y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=fish_odds$go_id, srt=45,adj=1, xpd=TRUE, cex=0.8)
	# y-axis
	# use automatic ticks for y-axis, but suppress below 0
	ti = par("yaxp") # min, max, n_intervals for automatic y-ticks
	ticks = seq(ti[1],ti[2],length.out=ti[3]+1)
	ticks = ticks[ticks>=0]
	axis(2, at=ticks, labels=10^ticks, cex.axis=0.8, las=2)
	# 95%-CI
	arrows(c(1:nrow(fish_odds),1:nrow(fish_odds)), log(c(fish_odds$CI_high,fish_odds$CI_low),10), c(1:nrow(fish_odds),1:nrow(fish_odds)), log(c(fish_odds$odds_ratio, fish_odds$odds_ratio),10), angle=90, code=1, length=0.03)
	# horizontal line at 1
	abline(h=0, col="red")
	axis(4,at=0,label="1", col="red", col.axis="red", las=1, cex.axis=0.8)
	# add significance asterix for fisher-test? - but does not account for multiple testing
	# rather FWER-over and under?

	# add pie charts
	space = limit/7
	radi_units = space/max(drawn_table[,5:6])
	for(i in 1:nrow(drawn_table)){
		b = drawn_table[i,5]
		w = drawn_table[i,6]
		add.pie(z=c(b,w), x=i, y=-limit/4, radius=(b+w)*radi_units, labels="")
	}
	text(x=1:nrow(fish_odds), y=par()$usr[3]+0.04*(par()$usr[4]-par()$usr[3]), labels=drawn_table[,5]+drawn_table[,6], xpd=TRUE, cex=0.8)
	
}


fisher = function(drawn, urn){
	# drawn white and black balls in node
	w = as.numeric(drawn["w"])
	b = as.numeric(drawn["b"])
	# white and black not drawn (left in urn)
	urn = urn_table[drawn["root_id"],]
	wn = urn["1"] - w
	bn = urn["0"] - b
	# contingency table
	conti = cbind(c(w,b),c(wn,bn))
	fish = fisher.test(conti)
	out = c(fish$estimate, fish$conf.int, fish$p.value) 
	return(out)	
}
