
# plot the odds-ratio and the 95%-CI from Fishers excact test (two-sided)
# show number of background and candidate genes in pie-charts with size dependent on total nr. of genes
# take table with candidate and background gene numbers prepared in plot_anno_scores

plot_hyper = function(aggrego, root_aggrego){
    
    # flip canidate and background columns
    colnames(aggrego)[c(2:3, 5:6)] = c("candi_genes", "bg_genes", "root_candi_genes", "root_bg_genes")
    colnames(root_aggrego)[c(2:3)] = c("candi_genes", "bg_genes")
    
    fish_odds = data.frame(aggrego, t(apply(aggrego[c(2:3, 5:6)], 1, fisher)))
    colnames(fish_odds)[9:12] = c("odds_ratio","ci95_low","ci95_high","p")

    # replace Inf or -Inf CI for plotting
    inf_go = fish_odds[is.infinite(fish_odds$odds_ratio),"go_id"]
    if(length(inf_go) > 0){
        warning(paste("The following GO-categories have an infinite odds ratio due to no background gene annotation: "), paste(inf_go, collapse=", "), sep="")
    }
    zero_go = fish_odds[fish_odds$odds_ratio==0,"go_id"]
    if(length(zero_go) > 0){
        warning(paste("The following GOs have an odds ratio of 0: "), paste(zero_go, collapse=", "), sep="")
    }

    ### plot the odds-ratios, with CI
    # TODO influenced by setting par(oma) before function call?
    op = par(no.readonly = TRUE) 
    # 3 panels
    layout(matrix(c(1,2,3,3),ncol=2), widths=c(4,1),heights=c(3,2))

    # odds-ratio
    par(mar=c(0.5,4,3,1), bty="l") #, bty="n")
    ymin = min(c(1,fish_odds$ci95_low[is.finite(fish_odds$ci95_low) & fish_odds$odds_ratio!=0]))
    ymax = max(c(1,fish_odds$ci95_high[is.finite(fish_odds$ci95_high) & fish_odds$odds_ratio!=0]))
    suppressWarnings(plot(fish_odds$odds_ratio, pch=19, ylab="", xaxt="n", xlab="", main="odds ratio", xlim=c(0.5,nrow(fish_odds)+0.5), ylim=c(ymin,ymax), cex.axis=0.8, log="y", las=2, 
    panel.first={grid(0, NULL, lty=1, col=colors()[2])}))
    # 95%-CI
    suppressWarnings(arrows(c(1:nrow(fish_odds),1:nrow(fish_odds)),c(fish_odds$ci95_high,fish_odds$ci95_low), c(1:nrow(fish_odds),1:nrow(fish_odds)), c(fish_odds$odds_ratio, fish_odds$odds_ratio), angle=90, code=1, length=0.03))
    # horizontal line at 1
    abline(h=1, col="#F15A60")
    axis(4,at=1,labels="1", col="#F15A60", col.axis="#F15A60", las=1, cex.axis=0.8)
    
    # pie charts
    par(mar=c(5.5,4,3,1), bty="l")
    plot(1,xlim=c(0.5,nrow(fish_odds)+0.5), ylim=c(0,1),type="n", main="annotated genes", xlab="", ylab="", xaxt="n", yaxt="n")
    radi_units = 0.4/log(max(rowSums(fish_odds[,2:3]))+1)
    for(i in 1:nrow(fish_odds)){
        w = fish_odds[i,2]
        b = fish_odds[i,3]
        add.pie(z=c(b,w), x=i, y=0.6, radius=log(b+w+1)*radi_units, labels="", col=c("#737373",fish_odds[i,"root_col"]))
    }
    text(x=1:nrow(fish_odds), y=0.08, labels=paste(fish_odds[,2],rowSums(fish_odds[,2:3]),sep=" / "),col= fish_odds$root_col, xpd=TRUE, cex=0.8)
    axis(1, at=1:nrow(fish_odds), labels=FALSE, cex.axis=0.8)
    text(x=1:nrow(fish_odds), y=-0.25, labels=fish_odds$go_id, srt=45, adj=1, xpd=TRUE, cex=0.8)
    
    # pie charts for root nodes
    par(mar=c(5.5,1,3,1), bty="l") # TODO, replace 3 with number of root nodes
    plot(1,xlim=c(-0.3,2), ylim=c(0.5,3.5),type="n", main="annotated genes\nroot nodes", xlab="", ylab="", xaxt="n", yaxt="n", yaxs="i")
    radi_units = 0.3/log(max(rowSums(root_aggrego[,2:3]))+1)
    for(i in 1:nrow(root_aggrego)){
        w = root_aggrego[i,2]
        b = root_aggrego[i,3]
        add.pie(z=c(b,w), x=1, y=i, radius=log(b+w+1)*radi_units, labels="", col=c("#737373",root_aggrego[i,"root_col"]))
    }
    text(x=1, y=(0.4 + 1:nrow(root_aggrego)), labels=root_aggrego$root_name, col=root_aggrego$root_col, cex=0.8)
    text(x=1, y=(-0.4 + 1:nrow(root_aggrego)), labels=paste(root_aggrego[,2],rowSums(root_aggrego[,2:3]), sep=" / "), cex=0.8, col=root_aggrego$root_col)
    
    par(op)
    
#    #TODO: also return p? add significance asterix for fisher-test? - but does not account for multiple testing, # rather FWER-over and under?
    # TODO: add go_name and root_name? (although this is easy with get_names)
    out = fish_odds[,-c(7,8)]
    return(invisible(out))
}


fisher = function(fish_line){
    # drawn white and black balls in node
    w = fish_line["candi_genes"]
    b = fish_line["bg_genes"]
    # white and black not drawn (left in urn)
    wn = fish_line["root_candi_genes"] - w
    bn = fish_line["root_bg_genes"] - b
    # contingency table
    conti = cbind(c(w,b),c(wn,bn))
    fish = fisher.test(conti)
    out = c(fish$estimate, fish$conf.int, fish$p.value)
    return(out)
}

