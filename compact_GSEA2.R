library(patchwork)
library(scales) # trans_new() is in the scales library

# function to transform in neglogscale the colorbar
neglog10_trans <- function(base = 10) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

compact_GSEA2 <- function(.GSEA, cutoff=0.05, title = "", truncate_label_at = 42) {
  .GSEA <- .GSEA@result %>%
    filter(qvalues <= cutoff) %>%
    arrange(NES) %>%
    select(Description, ID, NES, qvalues) %>%
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
    geom_rect(xmin=0, xmax=0, ymin=0,ymax=0, alpha = 0, aes(color = qvalues))+
    geom_col(aes(fill = NES > 0, alpha = qvalues), width = 0.75) +
    geom_text(aes(label = Description, x = x_pos, hjust = hjust), size = 2, color = "black") +
    geom_vline(xintercept = 0, linewidth = 0.25) +
    xlim(c(-max_NES, max_NES)) +
    scale_fill_manual(values = c("#3082BD", "#DE2C26")) +
    scale_alpha_continuous(range = c(0.5, 1), trans = "neglog10") +
    scale_color_gradient(low = "#8A8A8A80", high = "#8A8A8A", trans = "neglog10") + 
    theme_bw(8) +
    theme(panel.grid = element_blank(),
          # panel.border = element_rect(linewidth = 0.5, color="black", fill=NA),
          plot.margin = margin(t=0, b=0,unit="cm"),
          plot.title = element_text(size = 10, hjust = 0.5),
          plot.background = element_blank(),
          axis.line.y = element_blank(),
          panel.border = element_blank(),
          axis.line.x.bottom =  element_line(linewidth = 0.2, color="black"),
          axis.title.y = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 8, color = "black"),
          axis.ticks.length.y = unit(0, "null"),
          axis.text.y = element_blank(),
          text = element_text(size = 8, color = "black"),
          legend.position ="bottom",
          legend.text = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black"),
          legend.title = element_text(size = 8, color = "black", vjust = 0.96),
          legend.key.width = unit(1.5, "cm"),
          legend.key.height = unit(0.3, "cm"),
          #legend.title.align = 1,
          legend.direction = "horizontal"
          ) +
    xlab("Normalized Enrichment Score (NES)") +
    guides(fill = "none", alpha = "none",
           color = guide_colorbar(title.position = "left", ticks.colour =  "black", 
                                  frame.colour = "black",
                                  frame.linewidth = 0.2)) +
    labs(color = "adj. pvalue") + 
    ggtitle(title) +
    force_panelsizes(cols = unit(9, "cm"),
                     rows = unit((.GSEA %>% filter(qvalues <= cutoff ) %>% nrow())/2, "cm"))

plot(.plot)
# return a list, first item is the plot, second are the dimensions for saving it 
list("plot" = .plot,
     "height" = .GSEA %>% filter(qvalues <= cutoff) %>% nrow()/2 + 2)
}
