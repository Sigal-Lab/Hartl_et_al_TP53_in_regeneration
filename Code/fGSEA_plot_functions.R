# Functions taken and slightly adapted from fgsea package from CRAN, https://www.bioconductor.org/packages/release/bioc/html/fgsea.html
# Original copyright Alexey Sergushichev licensed under MIT+ license
myPlotEnrichment <- function (pathway, stats, gseaParam = 1, linewidth=1, group_pos = "Up", group_neg = "Down", color_pos = "red", color_neg = "blue", color_gene_bars = T)
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  
  so = sort(stats, decreasing=T)
  probs = c(1,0.95,0.9,0.75,0.5,0.25,0.1,0.05,0)
  percentiles = trunc(rank(so))/length(so)
  x_ticks = findInterval(-probs, -percentiles, all.inside=T) + 1
  x_ticks = sort(x_ticks)
  x_labels = prettyNum(so[x_ticks], digits=2)
  
  bottom_line <- min(bottoms)
  top_line <- max(tops)
  
  g1 <- ggplot(toPlot, aes(x = x, y = y)) + 
    geom_point(color = "green", size = 0.1) + 
    geom_hline(yintercept = top_line, colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = bottom_line, colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = 0, colour = "black") + 
    geom_line(color = "green", linewidth=linewidth) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle = 30, hjust = 1)) + #,
    #axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks = element_blank() ) +
    coord_cartesian(xlim=c(0,max(toPlot$x))) + 
    scale_x_continuous(breaks=x_ticks, labels=x_labels, sec.axis = sec_axis(~.+10)) + 
    scale_y_continuous(breaks=seq(0,max(tops), 0.2)) 
  
  if(color_gene_bars) {
    g1 <- g1 + geom_segment(data = data.frame(x = pathway, v = so[pathway]), 
                 mapping = aes(x = x, 
                               y = bottom_line - 2.2*diff/2, 
                               xend = x, 
                               yend = bottom_line - 1.1*diff/2,
                               color = v), 
                 size = 0.2) +
      scale_color_gradient2(low="blue", high="red", mid="grey", midpoint = 0)
  } else {
    g1 <- g1 + geom_segment(data = data.frame(x = pathway, v = so[pathway]), 
                 mapping = aes(x = x, 
                               y = bottom_line - 2.2*diff/2, 
                               xend = x, 
                               yend = bottom_line - 1.1*diff/2), 
                 size = 0.2)
  }

    #geom_tile(data=data.frame(x = pathway), stat = "identity", aes(x=x,y=1, fill=x))  + scale_fill_gradient2(low="blue", high="red", mid="grey") +
  g1 <- g1 + labs(x = "Ranks based on moderated t statistic", y = "Enrichment score") + 
    annotate("text", x=0, y = bottom_line - (diff*0.1), label= paste0('paste(bold(',group_pos,'))'), color=color_pos, parse=T, hjust="inward", vjust="top") + 
    annotate("text", x=length(so), y = bottom_line - (diff*0.1), label= paste0('paste(bold(',group_neg,'))'), color=color_neg, parse=T, hjust="inward", vjust="top") +
    theme(legend.position="none")
  
  g1
}
