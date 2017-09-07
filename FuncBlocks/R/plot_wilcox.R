
# plot the distribution of annotated scores for nodes and root nodes
# take table with sums of scores per node from plot_anno_scores
# create output with list of 3 tables(node-anno, root-anno, node->root-mapping)

plot_wilcox = function(anno_scores, root_anno_scores, root_info){
	
	colnames(anno_scores)[3] = "score"
	colnames(root_anno_scores)[3] = "score"
	
	# go_id -> root table
	root_names = get_names(anno_scores[,1])[,3] # get names
	node_to_root = data.frame(go_id=anno_scores[,1], root_id=root_info[match(root_names, root_info[,3]), 1])

	# layout
	layout(matrix(c(1,2),ncol=2), widths=c(5,2))
	ylim = range(c(anno_scores$score, root_anno_scores$score))
#	ylim = range(c(anno_scores$score, root_anno_scores$score)/ mean(root_cols$median)) # TODO: proper range if normalized
#	yrange = ylim[2] - ylim[1]
#	ylim[2] = ylim[2] + 0.3*yrange
	op = par(no.readonly = TRUE) 
	# plot GO-categories
	par(mar=c(6.5,4,4,2), bty="l") #, bty="n") # mar default=c(5, 4, 4, 2)
	violin(anno_scores, root_info, node_to_root, ylim)
	# plot root nodes
	par(mar=c(6.5,1,4,2), bty="l")
	violin(root_anno_scores, root_info, node_to_root, ylim, root=TRUE)
	par(op) 

	out = list(node_anno=anno_scores, root_anno=root_anno_scores, node_to_root=node_to_root)

	return(invisible(out))
}


## plot multiple violins side by side
violin = function(plotty, root_info, node_to_root, ylim, root=FALSE){
	# find root-node for every GO and add median and color (gos=[node, root] for every node to plot)
	if (root){
		gos = data.frame(unique(plotty[,1]), unique(plotty[,1]))
		mainy = "root nodes"
	} else {
		gos = node_to_root[match(unique(plotty[,1]), node_to_root[,1]) ,]
		mainy = "scores of annotated genes"
	}
	# root_info: [id, median, name, col] of root for every node
	root_info = root_info[match(gos[,2], root_info[,1]), ] 
	# xlim
	xoff = 0.5
	xlim = c(1-xoff, nrow(gos)+xoff)
	# plot
	plot(1, type="n", ylim=ylim, xlim=xlim, ylab="score", xaxt="n", xlab="", main=mainy) #, log="y")
	if(root){
		xlabel = root_info$root_name
		labelcol = root_info$root_col
	} else {
		xlabel = gos[,1]
		labelcol = "black"
		# plot horizontal line at root median
		abline(h=root_info[,2], col=root_info$root_col)
	}
	# x-axis
	xlabpos = ylim[1]-(ylim[2]-ylim[1])/12
#	xlabpos = ylim[1]-(ylim[2]-ylim[1])/50  # for log-axis
	axis(1, at=1:nrow(gos), labels=FALSE, cex.axis=0.9)
	text(x=1:nrow(gos), y=xlabpos , labels=xlabel, srt=45, adj=c(1,1), xpd=TRUE, cex=0.8, col=labelcol)
	for(i in 1:nrow(gos)){
		scores = plotty[plotty$go_id == gos[i,1], "score"]
		# divide by root median
#		scores = scores / gos[i,"median"]
		if (length(scores)==1) {
			points(x=i, y=scores, col=root_info[i,"root_col"], pch=16, cex=1.5)
		} else {
			vioplot(scores, at=i, col=root_info[i,"col"], add=TRUE)
		}
#		points(i+runif(length(scores),-0.15, 0.15), scores)
	}
	# n genes
	n = table(plotty$go_id)
	text(x=1:nrow(gos), y=ylim[2], labels=paste0("n=",n[gos[,1]]), xpd=TRUE, cex=0.8, col=labelcol, pos=3, offset=0.5)
}
