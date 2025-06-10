library(enrichplot)
library(glue)
library(scales)
options(scipen = -1) # turn on scientific notation after 3 decimals

multiPathwayGSEAplot <- function(GSEAresult, geneSetID, genesAlpha = 0.5,
                           simplifyCurve = TRUE, linesColors = c("#297FB9", "#C03A2B", "#25AF60", "#8F44AD", "#2D3E50",  "#F49C12")) {

    custom_theme <- function() {
    #theme_void(8) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.33, lineend="butt"),
          axis.ticks.length=unit(.07, "cm"),
          axis.line = element_line(color = "black", linewidth = 0.33, lineend="square"),
          axis.text = element_text(size = 8, color = "black"),
          axis.title.x = element_text(color = "black", size = 10),
          axis.title.y = element_text(color = "black", size = 10, angle = 90),
          plot.margin = margin(t=0, b=0,unit="cm"),
          axis.ticks.margin = unit(0, "null"),
          text = element_text(size = 8),
          legend.position = "top",
          legend.justification = c(0, 0.5),
          legend.margin = margin(0, 0, 0, 0),   
          legend.box.just = "left",
          legend.text  = element_text(size = 8),
          legend.background = element_blank(),
          legend.key = element_blank()) 
  }
  
  # extract gsea plotting raw data for the signatures of interest
  extract_gsdata <- function(.GSEAresult, .geneSetID) {
    geneSetID.in_result <- .geneSetID[.geneSetID %in% names(.GSEAresult@geneSets)]
    gsdata <- map(geneSetID.in_result, ~enrichplot:::gsInfo(.GSEAresult, .)) %>%
      bind_rows() 
  }
  
  # handle signatures coming from different GSEA analysis
  if(length(GSEAresult) == 1) {
    gsdata <- extract_gsdata(unlist(GSEAresult), geneSetID)
  } else {
    gsdata <- map(GSEAresult, ~extract_gsdata(.x, geneSetID)) %>%
      bind_rows()
  }

  if(simplifyCurve){
    gsdata <- gsdata %>%
      group_by(Description) %>%
      # remove points by gsdata for which the previous point has
      # an higher runningScore value and the next point has a lover runningScore value
      filter(!(runningScore < dplyr::lag(runningScore) & runningScore > dplyr::lead(runningScore)) |
               is.na(dplyr::lag(runningScore)) | is.na(dplyr::lead(runningScore))) %>%
      ungroup()
  }
  
  # Variables -------
  #enrichmentScore <- GSEAresult@result[geneSetID, "enrichmentScore"]
  y_axis.range = max(gsdata$runningScore) - min(gsdata$runningScore)
  y_axis.min <- min(gsdata$runningScore) - 0.05*y_axis.range
  y_axis.max <-  max(gsdata$runningScore) + 0.05*y_axis.range
  
  gsdata <- gsdata %>%
    # compute y position for each enrich line
    mutate(Description = factor(Description, levels = geneSetID) %>% fct_drop(),
    # yposition adds a margin for x axis title and lables + additional margin for each pavlue label
           yend = y_axis.min - 0.35*y_axis.range - 0.3*as.numeric(Description)*y_axis.range,
           y =  yend - 0.047*y_axis.range)
  
  # compute gene set p-value statistic --------------------------
  label <- GSEAresult %>%
    # in case of multiple GSEA results, combine them
    {if(typeof(.) == "list") {map(., as.data.frame) %>% bind_rows()} else {as.data.frame(.)}} %>%
  #  map(as.data.frame) %>% bind_rows() %>%
    filter(Description %in% geneSetID) %>%
    rename_with(~str_replace(., "qvalue", "qvalues"), any_of("qvalue")) %>%
    mutate(label = glue("NES: {NES}; p: {pvalue}; padj: {qvalue}",
                        NES = signif(NES, 2),
                        pvalue = signif(pvalue, 2),
                        qvalue = signif(qvalues, 2)
    )) %>%
    # add y position
    mutate(Description = factor(Description, levels = geneSetID) %>% fct_drop(),
          # yposition adds a margin for x axis title and lables + additional margin for each pavlue label
          y = y_axis.min - 0.23*y_axis.range - 0.3*as.numeric(Description)*y_axis.range)
  
  # Extract points for max Enrichment score
  ES.data <- gsdata %>%
    group_by(Description) %>%
    slice_max(order_by = abs(runningScore), n=1, with_ties = F) %>%
    select(x, Description, runningScore) %>%
    ungroup()
  
  
  # Plot ----------
  ggplot() + 
    coord_cartesian(clip = "off", ylim = c(y_axis.min, y_axis.max), expand = 0) +
    # add zero line
    geom_hline(yintercept = 0, color = "#a9aaaa", linetype = "dashed", linewidth = 0.33) +
    # plot max ES sline
    geom_segment(data = ES.data, aes(x = x, xend = x, y = 0, yend = runningScore, color = Description),
                 linetype = "dashed", linewidth = 0.33) +
    # plot enrichment score line
    geom_line(data = gsdata, aes(x=x, y = runningScore, color = Description), linewidth = 2/.pt) +
    # plot genes in the set
    geom_segment(data = gsdata %>% filter(position == 1),
                 aes(x = x, xend=x, yend = yend, y =  y, color = Description),
                 size = 0.3, alpha = genesAlpha) +
    # add statistics
    geom_text(data =label, aes(x = 0, y = y, label = label), inherit.aes = F,
              size = 10/.pt, hjust = 0) +
    custom_theme() +
    xlab("Gene Rank") + ylab("Enrichment Score") +
    scale_y_continuous(labels = label_number(drop0trailing = TRUE)) +
    scale_x_continuous(breaks = seq(0, 55000, by = 5000), limits = c(0, max(gsdata$x))) +
    scale_color_manual(values = linesColors) +
    labs(color = NULL) + # remove legend title
    guides(color = guide_legend(keywidth = unit(0.3, "cm"), ncol = 2)) + # shorten legend line symbol
    plot_layout(heights = unit(2.13, "cm"), widths = unit(3, "cm"))
}
