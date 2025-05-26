library(scales)

custom_gseaplot2 <- function(x, geneSetID, genesAlpha = 0.5, simplifyCurve = TRUE) {
  statistics <- x %>% as.data.frame() %>% filter(Description == geneSetID)
  
  gsdata <- enrichplot:::gsInfo(x, geneSetID)
  
  if(simplifyCurve){
    gsdata <- gsdata %>%
      # remove points by gsdata for which the previous point has
      # an higher runningScore value and the next point has a lover runningScore value
      filter(!(runningScore < dplyr::lag(runningScore) & runningScore > dplyr::lead(runningScore)) |
               is.na(dplyr::lag(runningScore)) | is.na(dplyr::lead(runningScore)))
  }
  
  # statistics label
  signif_label <- ifelse(nrow(x@result) > 1,
                         paste0("NES: ", round(statistics$NES,2),"\np: ", format(statistics$pvalue, digits=2, scipen = 0), "\nqval: ",  format(statistics$qvalues, digits=2, scipen = 0)),
                         paste0("NES: ", round(statistics$NES,2),"\np: ", format(statistics$pvalue, digits=2, scipen = 0))
  )
  
  
  enrichmentScore <- x@result[geneSetID, "enrichmentScore"]
  
  # identify the datapoint that leads to the ES
  es.df <- gsdata %>%
    slice_max(order_by = abs(runningScore), n=1, with_ties = F)
  
  # verify its matching the enrichmentScore
  if(abs(es.df$runningScore - enrichmentScore) > 0.1) {
    stop("Enrichment Score does not match extracted runningScore")
  }
  
  y_axis.range = max(gsdata$runningScore) - min(gsdata$runningScore)
  y_axis.min <- min(gsdata$runningScore) - 0.05*y_axis.range
  y_axis.max <-  max(gsdata$runningScore) + 0.05*y_axis.range
  
  # setup theme
  p1 <- ggplot() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.33, lineend="butt"),
          axis.ticks.length=unit(.07, "cm"),
          axis.line = element_line(color = "black", linewidth = 0.33, lineend="square"),
          axis.text = element_text(size = 8, color = "black"),
          axis.title = element_text(color = "black", size = 10),
          axis.title.x = element_text(color = "black", size = 10, margin = margin(t=0.2, b=0,unit="cm")),
          plot.margin = margin(t=0, b=0,unit="cm"),
          text = element_text(size = 8),
          legend.position = "none") +
    scale_y_continuous(labels = label_number(drop0trailing = TRUE)) +
    scale_x_continuous(breaks = seq(0, 55000, by = 5000), limits = c(0, max(gsdata$x))) +
    #scale_x_continuous(expand=c(0,0)) +
    coord_cartesian(clip = "off", ylim = c(y_axis.min, y_axis.max), expand = 0) +
    # plot genes in the set
    geom_segment(data = gsdata %>% filter(position == 1), aes(x = x, xend=x, yend = y_axis.min - 0.2*y_axis.range, y =  y_axis.min - 0.247*y_axis.range, ymax =), color="black", linewidth = 0.3, alpha = genesAlpha) +
    # plot enrichment score line
    geom_segment(data=es.df, aes(x = x, xend = x, y = 0, yend=runningScore),
                 colour = "#2479ae", linetype = "dashed", linewidth = 0.33) +
    # geom_segment(data=es.df, aes_(x = -Inf, xend = ~es, y = enrichmentScore, yend=enrichmentScore),
    #              colour = "#2479ae", linetype = "dashed", linewidth = 0.33) +
    geom_hline(yintercept = 0, color = "#a9aaaa", linetype = "dashed", linewidth = 0.33) +
    # add area
    # geom_polygon(aes(y = runningScore),
    #              size=1, fill="#1A9F62", alpha = 0.3) +
    # add line
    geom_line(data = gsdata, aes(x=x, y = runningScore),
              size=2/.pt, color = "#2479ae") +
    # position for cohordinates based on positive or negative enrichment score
    {if(enrichmentScore > 0)annotate("text",
                                     x=max(gsdata$x)*0.97, y=enrichmentScore*0.97,
                                     hjust = "right", vjust = "top",size=3,
                                     label = signif_label)} +
    {if(enrichmentScore < 0)annotate("text",
                                     x=max(gsdata$x)*0.03, y=enrichmentScore*0.97,
                                     hjust = "left", vjust = "bottom",size=3,
                                     label = signif_label)} +
    ylab("Enrichment Score") + xlab("Gene Rank") +
    ggtitle(geneSetID)
  
  p1 + plot_layout(heights = unit(2.13, "cm"), widths = unit(3, "cm")) 
  #ggh4x::force_panelsizes(rows = unit(2.13, "cm"), cols = unit(3, "cm"), respect = F)
}
