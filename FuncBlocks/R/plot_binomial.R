
# plot the p(A/B) and the 95%-CI from binomial test
# take table with sums of scores per node from plot_anno_scores


plot_binomial = function(aggrego, root_aggrego){
	
	# p_A for root nodes
	colnames(aggrego)[c(2,3,5,6)] = c("A_node", "B_node", "A_root", "B_root")
	aggrego$p_A_root = aggrego[,5] / (aggrego[,5] + aggrego[,6])
	root_aggrego$p_A = root_aggrego[,2] / (root_aggrego[,2] + root_aggrego[,3])
	# perform binomial test for every GO
	binom = data.frame(aggrego, t(apply(aggrego[,c(2,3,9)], 1, binomy)))
	colnames(binom)[c(10:13)] = c("p_A_node","ci95_low","ci95_high","p")
	
	### plot p_A_node, with CI
	op = par(no.readonly = TRUE) 
	# 3 panels
	layout(matrix(c(1,2,3,3),ncol=2), widths=c(4,1),heights=c(3,2))

	# p_A_node
	par(mar=c(0.5,4,3,1), bty="l") #, bty="n")
	rangy = range(c(binom$ci95_low, binom$ci95_high, binom$p_A_root))
	ymin = max(0, rangy[1] - 0.2 * (rangy[2] - rangy[1]))
	ymax = min(1, rangy[2] + 0.2 * (rangy[2] - rangy[1]))
#	ymin = 0
#	ymax = 1
	suppressWarnings(plot(binom$p_A_node, pch=19, ylab="", xaxt="n", xlab="", main="p(A) in node", xlim=c(0.5,nrow(binom)+0.5), ylim=c(ymin,ymax), cex.axis=0.8, las=2, 
	panel.first={grid(0, NULL, lty=1, col=colors()[2])}))
	# 95%-CI
	suppressWarnings(arrows(c(1:nrow(binom),1:nrow(binom)),c(binom$ci95_high,binom$ci95_low), c(1:nrow(binom),1:nrow(binom)), c(binom$p_A_node, binom$p_A_node), angle=90, code=1, length=0.03))
	# horizontal line for root-node p_A
	abline(h=root_aggrego$p_A, col=root_aggrego$root_col)
	
	## pie charts
	par(mar=c(5.5,4,3,1), bty="l")

	plot(1,xlim=c(0.5,nrow(binom)+0.5), ylim=c(0,1),type="n", main="proportion of A", xlab="", ylab="", xaxt="n", yaxt="n")
	radi_units = 0.4/log(max(rowSums(binom[,2:3]))+1)
	for(i in 1:nrow(binom)){
		a = binom[i,2]
		b = binom[i,3]
		add.pie(z=c(b,a), x=i, y=0.6, radius=log(a+b+1)*radi_units, labels="", col=c("#737373", binom[i,"root_col"]))
	}
	text(x=1:nrow(binom), y=0.08, labels=paste(binom[,2], rowSums(binom[,2:3]),sep=" / "), col=binom$root_col, xpd=TRUE, cex=0.8)
	axis(1, at=1:nrow(binom), labels=FALSE, cex.axis=0.8)
	text(x=1:nrow(binom), y=-0.25, labels=binom$go_id, srt=45, adj=1, xpd=TRUE, cex=0.8)
	
	# pie charts for root nodes
	par(mar=c(5.5,1,3,1), bty="l")
	plot(1,xlim=c(-0.3,2), ylim=c(0.5,3.5),type="n", main="proportion of A\nroot nodes", xlab="", ylab="", xaxt="n", yaxt="n", yaxs="i")
	radi_units = 0.3/log(max(rowSums(root_aggrego[,2:3]))+1)
	for(i in 1:nrow(root_aggrego)){
		a = root_aggrego[i,2]
		b = root_aggrego[i,3]
		add.pie(z=c(b,a), x=1, y=i, radius=log(b+a+1)*radi_units, labels="", col=c("#737373",root_aggrego[i,"root_col"]))
	}
	text(x=1, y=(0.4 + 1:nrow(root_aggrego)), labels=root_aggrego$root_name, col=root_aggrego$root_col, cex=0.8)
	text(x=1, y=(-0.4 + 1:nrow(root_aggrego)), labels=paste(root_aggrego[,2],rowSums(root_aggrego[,2:3]), sep=" / "), cex=0.8, col=root_aggrego$root_col)
	
	par(op)

	out = binom[,-c(7,8)]

	return(invisible(out))
}


binomy = function(bino_line){
	# binom.test(c(success, failure), p)
	bino_line = unname(bino_line)
	bino = binom.test(c(bino_line[1], bino_line[2]), p=bino_line[3])
	out = c(bino$estimate, bino$conf.int, bino$p.value) 
	return(out)	
}

