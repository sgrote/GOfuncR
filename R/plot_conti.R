
# plot the odds-ratio and the 95%-CI from Fishers excact test (two-sided)
# show A,B,C,D counts in pie-charts with size dependent on total nr. of genes
# take table with sums of scores per node from plot_anno_scores

# TODO: test plotted odds ratio with an example

plot_conti = function(aggrego){
    
    # don't need root node for that
    aggrego = aggrego[,1:5] # TODO: maybe avoid adding root-info in plot_anno_scores in the first place?
    
    ### perform fishers exact test for every GO
    fish_odds = data.frame(aggrego, t(apply(aggrego, 1, fisher_conti)))
    colnames(fish_odds)[c(2:9)] = c("A","B","C","D","odds_ratio","ci95_low","ci95_high","p")
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
#   pie_cols = colors()[c(592,566,503,556)]

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

