library(patchwork)
library(scales) # trans_new() is in the scales library
library(ggh4x)


# function to transform in neglogscale the colorbar
neglog10_trans <- function(base = 10) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(name = paste0("reverselog-", format(base)),
            transform = trans,
            inverse = inv, 
            breaks = log_breaks(base = base, n = 5), 
            domain = c(1e-100, Inf))
}

text_size = 10
title_size = 11
linewidth = 0.7/2.141959


gseaBarplot <- function(.GSEA, cutoff=0.05, title = "", truncate_label_at = 40) {
  .GSEA <- .GSEA@result %>%
    filter(qvalue <= cutoff) %>%
    arrange(NES) %>%
    select(Description, ID, NES, qvalue) %>%
    mutate(Description = str_trunc(Description, width = truncate_label_at, side = "right", ellipsis = "..")) %>%
    # specify which side plot the labels
    mutate(hjust = ifelse(NES > 0, 1, 0),
           # add some distance form x axis
           x_pos = ifelse(NES > 0, -0.05*max(abs(.$NES)), 0.05*max(abs(.$NES))))
  
  .GSEA <- mutate(.GSEA) %>%
    mutate(ID = factor(ID, levels = c(.GSEA$ID)))
  
  max_NES <- max(abs(.GSEA$NES))
  
  .plot <-  ggplot(.GSEA, aes(x=NES, y=ID)) +
    # mock color gradient to supplement continous alpha scale not available in ggplot
    geom_rect(xmin=0, xmax=0, ymin=0,ymax=0, alpha = 0, aes(color = qvalue))+
    geom_col(aes(fill = NES > 0, alpha = qvalue), width = 0.75) +
    geom_text(aes(label = Description, x = x_pos, hjust = hjust), size = 8/.pt, color = "black") +
    geom_vline(xintercept = 0, linewidth = linewidth) +
    xlim(c(-max_NES, max_NES)) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("FALSE" ="#2479Ae", "TRUE" = "#C93E52")) +
    scale_alpha_continuous(range = c(0.5, 1), trans = "neglog10") +
    scale_color_gradient(low = "#8A8A8A80", high = "#8A8A8A", trans = "neglog10") + 
    theme_bw(8) +
    theme(panel.grid = element_blank(),
          # panel.border = element_rect(linewidth = 0.5, color="black", fill=NA),
          plot.margin = margin(t=0, b=0,unit="cm"),
          plot.title = element_text(size = 11, hjust = 0.5),
          plot.background = element_blank(),
          axis.line.y = element_blank(),
          panel.border = element_blank(),
          axis.line.x.bottom =  element_line(linewidth = linewidth, color="black"),
          axis.title.y = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = linewidth),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 11, color = "black"),
          axis.ticks.length.y = unit(0, "null"),
          axis.text.y = element_blank(),
          text = element_text(size = 8, color = "black"),
          legend.position ="bottom",
          legend.text = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, color = "black"),
          legend.title = element_text(size = 10, color = "black", vjust = 0.96),
          legend.key.width = unit(0.75, "cm"),
          legend.key.height = unit(0.3, "cm"),
          #legend.title.align = 1,
          legend.direction = "horizontal"
    ) +
    xlab("Normalized Enrichment Score (NES)") +
    guides(fill = "none", alpha = "none",
           color = guide_colorbar(title.position = "left", ticks.colour =  "black", 
                                  frame.colour = "black",
                                  frame.linewidth = linewidth)) +
    labs(color = "adj. pvalue") + 
    ggtitle(title) +
    force_panelsizes(cols = unit(2.857*2, "cm"),
                     rows = unit((.GSEA %>% filter(qvalue <= cutoff ) %>% nrow())/2.5, "cm"))
  
  plot(.plot)
  # return a list, first item is the plot, second are the dimensions for saving it 
  list("plot" = .plot,
       "height" = .GSEA %>% filter(qvalue <= cutoff) %>% nrow()/2.5 + 2)
}

