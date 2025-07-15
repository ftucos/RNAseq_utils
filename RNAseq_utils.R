library(msigdbr)
library(clusterProfiler)
library(patchwork)
library(ggrepel)
library(egg)
library(gridExtra)
library(tidyverse)
library(colorspace)
library(scales)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")

## px to pt conversion factor
.px = .pt*72.27522/96

# Annotate Result Table ---------------------------
ResToTable <- function(res, package) {
  if (package == "DESeq2") {
    # res obtained from DESeq2::result(dds) or DESeq2::lfcShrink(dds) functions 
    result.table <- res %>% 
      as.data.frame() %>%
      rownames_to_column("ensembl_gene_id") %>%
      left_join(ensembl2symbol) %>%
      select(external_gene_name, log2FoldChange, pvalue, padj, mean_gene_abundance = baseMean, gene_biotype, description, ensembl_gene_id)
    
  } else if (package == "edgeR") {
    # res obtained from edgeR::topTags(qlf) function
    result.table <- res$table %>%
      rownames_to_column("ensembl_gene_id") %>%
      left_join(ensembl2symbol) %>%
      select(external_gene_name, log2FoldChange = logFC, pvalue = PValue,  padj = FDR, mean_gene_abundance = logCPM, gene_biotype, description, ensembl_gene_id)
  } else {
    warning("package argument must be either 'DESeq2' or 'edgeR'!")
    return(NULL)
  }
}

# Volcano Plot -------------------------------------
volcanoPlot <- function(result, title = element_blank(), thrLog2FC = 1, thrPadj = 0.05, pctLabel = 0.005, protein_coding_label_only = FALSE) {
  
  # requires you to have applied ResToTable()
  x_limits <- c(-max(abs(result$log2FoldChange), na.rm=T), max(abs(result$log2FoldChange), na.rm=T))
  
  label_cutof <- result %>%
    # plot labels only for genes within the threshold defined
    filter(padj <= thrPadj, abs(log2FoldChange) >= thrLog2FC) %>%
    # rank priority of genes to plot based on combination of log2FC and padj
    # annotate protein coding only 
    filter((!protein_coding_label_only) | (ensembl_gene_id %in% (ensembl2symbol %>% filter(gene_biotype == "protein_coding") %>% pull("ensembl_gene_id")))) %>%
    mutate(DEscore = abs(log2FoldChange * -log10(pvalue))) %>%
    pull(DEscore) %>%
    # extract only the fraction of labels to plot
    quantile(1-pctLabel)
  
  df <- result %>%
    # Color code base on up/downregulation state defined by custom thresholds
    mutate(color = case_when(
      log2FoldChange >= thrLog2FC & padj <= thrPadj ~ "up",
      log2FoldChange <= -thrLog2FC & padj <= thrPadj ~ "down"),
      # labels to be shown need to be matching upregulated genes and be within the filtering threshold
      lab = ifelse(
        # log2 Fold Change is higher than selected threshold
        abs(log2FoldChange) >= thrLog2FC &
          # padj is below selected threshold  
          padj <= thrPadj &
          # annotate protein coding only 
          ((!protein_coding_label_only) | gene_biotype == "protein_coding") &
          abs(log2FoldChange * -log10(pvalue)) >= label_cutof,
        external_gene_name, NA)
    ) %>%
    # plot on background genes that are not within the filtering threshold
    arrange(color)
  
  
  
  ggplot(df, aes(x=log2FoldChange, y=-log10(padj), color = color)) +
    geom_point(size = 1) +
    geom_vline(xintercept = c(-thrLog2FC, thrLog2FC), linetype = "dashed")+
    geom_hline(yintercept = -log10(thrPadj), linetype = "dashed")+
    geom_text_repel(data = df %>% filter(!is.na(lab)),
                    aes(label = lab), size=3, color="black", segment.alpha = 0.7, max.overlaps = 40)+
    theme_bw(base_line_size = 0.75/.px, base_size = 10) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 2*0.75/.px),
          legend.position = "none") +
    xlim(x_limits) +
    scale_color_manual(values = c("dodgerblue", "brown1"), na.value = "grey") +
    xlab(expression(log["2"](Fold~Change))) + ylab(expression(-log["10"](Adjusted~p~Value))) + 
    ggtitle(title)
}

## Prepare imput data for GSEA
PrepareForGSEA <- function(result, rankingMetric = c("log2FoldChange", "signed_pvalue", "combined_score"), .na.rm = FALSE) {
  
  
  # this will (a) pick the first entry of the vector as the default
  # (b) partial‐match or reject invalid values
  rankingMetric <- match.arg(rankingMetric)
  
  if (rankingMetric == "log2FoldChange") {
    # chose whether remove genes with log2FC undefined (lowly expressed) or set theire log2FC to 0 since they are not activated by the treatment
    if (.na.rm == TRUE) {
      result <- result %>%
        filter(!is.na(log2FoldChange))
    } else {
      result <- result %>%
        mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange))
    }
    
    # Generate a sorted named vector with ensembl as name and log2FC as value
    .input <- setNames(result$log2FoldChange, result$ensembl_gene_id) %>%
      sort(decreasing = TRUE)
    
  } else if (rankingMetric == "signed_pvalue") {
    if (.na.rm == TRUE) {
      result <- result %>%
        filter(!is.na(pvalue), !is.na(log2FoldChange))
    } else {
      result <- result %>%
        mutate(pvalue = ifelse(is.na(pvalue), 1, pvalue),
               log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange))
    }
    
    # Generate a sorted named vector with ensembl as name and log2FC as value
    .input <- setNames(sign(result$log2FoldChange)*-log10(result$pvalue), result$ensembl_gene_id) %>%
      sort(decreasing = TRUE)
    
  } else if (rankingMetric == "combined_score") {
    result <- result %>%
      mutate(combined_score = log2FoldChange * -log10(pvalue))
    
    if (.na.rm == TRUE) {
      result <- result %>%
        filter(!is.na(combined_score))
    } else {
      result <- result %>%
        mutate(combined_score = ifelse(is.na(combined_score), 0, combined_score))
    }
    
    # Generate a sorted named vector with ensembl as name and log2FC as value
    .input <- setNames(result$combined_score, result$ensembl_gene_id) %>%
      sort(decreasing = TRUE)
    
  } else {
    stop("Invalid ranking metric. Choose either 'log2FoldChange', 'pvalue', or 'combined_score'.")
  }
  
  return(.input)
}

## Plot GSEA table
plotGSEA <- function(.GSEA, title = "", cutoff = 0.05, subgroup = "all", fixed_dimensions=F) {
  
  if (subgroup == "pos") {
    .GSEA <- .GSEA %>% filter(NES > 0)
  } else if (subgroup == "neg") {
    .GSEA <- .GSEA %>% filter(NES < 0)
  }
  plot <- .GSEA %>%
    filter(qvalues <= cutoff) %>%
    mutate(Description = factor(Description, levels = 
                                  .GSEA %>% arrange(NES) %>% pull(Description))) %>%
    ggplot(aes(x=Description, y=NES, fill=qvalues)) +
    geom_col() + 
    theme_bw(base_line_size = 0.75/.px, base_size = 10) +
    scale_fill_continuous(type = "gradient",
                          low = "#E74C3C", high = "#F1C40F",
                          space = "Lab", na.value = "grey70", guide = "colourbar")  +
    coord_flip() +
    ggtitle(paste0("GSEA: ", title)) +
    xlab("") +
    ylab("Normalized Enrichment Score") +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot")
  
  # Fixed dimensions based on number of pathways
  if(fixed_dimensions) {
    plot <-  set_panel_size(plot,
                            # define plot hieght based on the numer of pathways to be plotted
                            height=unit((.GSEA %>% filter(qvalues <= cutoff ) %>% nrow())/1.5, "cm"),
                            width=unit(8, "cm")
    )
    # plot the plot
    plot(plot)
    # return the plot (for export preservig dimensions)
    plot
  } else {
    # just plot the plot
    plot(plot)
    # and return it for export
    plot
  }
}

## Run GSEA
runGSEA <- function(result, annotation, rankingMetric = "log2FoldChange", title = "", cutoff = 0.05, plot = FALSE, toDataFrame = TRUE, custom_annotation=F, ...) {
  if(custom_annotation) {
    # custom t2g object to be passed into annotatation argument
    selected_t2g <- annotation
  } else {
    selected_t2g <- # extract t2g of interest based on string passed into annotation
      t2g %>% filter(gs_subcat == annotation) %>% select(gs_name, ensembl_gene)
  }
  
  .em <- GSEA(PrepareForGSEA(result, rankingMetric),
              TERM2GENE = selected_t2g,
              pAdjustMethod = "fdr",
              by = "fgsea",
              minGSSize = 5, maxGSSize = 2000,
              pvalueCutoff = 1,
              nPermSimple = 10000,
              #nPerm = 10000,
              eps = 0,
              seed = FALSE # use session seed
  )
  .em.summary <- as.data.frame(.em)
  
  # If not specified, use the object name as plot title
  title <- ifelse(title == "", deparse(substitute(result)), title)
  # PLot
  if (plot == TRUE) {
    plot(plotGSEA(.em.summary, title = title, cutoff = cutoff, ...))
  }
  # return the summary table
  if(toDataFrame == TRUE) {
    return(.em.summary)
  } else {
    return(.em)
  }
}


# ORA -----------------------------------
## Plot ORA table
plotORA <- function(.ORA, title = "", cutoff = 0.05, subgroup = "all", truncate_label_at = 45) {
  .ORA <- .ORA %>% 
    filter(qvalue <= cutoff) %>%
    # Truncate long strings
    mutate(Description = str_trunc(Description, width = truncate_label_at, side = "right", ellipsis = ".."))
  
  # use it to set a common range scale
  max_enrichment <- .ORA %>% filter(qvalue <= cutoff) %>% pull("EnrichmentRatio") %>% max()
  
  .ORA.pos = .ORA %>% filter(direction == "Up")
  .ORA.neg = .ORA %>% filter(direction == "Down")
  
  plot.pos <- .ORA.pos %>%
    mutate(neglog10p = -log10(qvalue),
           Description = factor(Description, levels = 
                                  .ORA.pos %>% arrange(EnrichmentRatio) %>% pull(Description) %>% unique())) %>%
    ggplot(aes(x=Description, y=EnrichmentRatio, fill=neglog10p, label=Description)) +
    geom_col() + 
    geom_text(y=max_enrichment*0.05, hjust = 0, size = 3.5)+
    scale_y_continuous(limits = c(0, max_enrichment))+
    theme_bw(base_line_size = 0.75/.px, base_size = 10) +
    scale_fill_continuous(type = "gradient",
                          high = "#DC3D2D", low = "#FED98B",
                          space = "Lab", na.value = "grey70",
                          label =  ~round(., digits = 1),
                          guide = "colourbar", name = "-Log10(adj p-value)")  +
    coord_flip() +
    ggtitle("ORA: enrichment of upregulated genes") +
    xlab("") +
    ylab("Enrichment Ratio") +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black", linewidth = 0.75/.px),
          axis.text = element_text(color = "black", size = 10),
          legend.text = element_text(color = "black", size = 10),
          plot.title.position = "plot",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          title = element_text(size = 10)
          )
  
  plot.neg <- .ORA.neg %>%
    mutate(neglog10p = -log10(qvalue),
           Description = factor(Description, levels = 
                                  .ORA.neg %>% arrange(EnrichmentRatio) %>% pull(Description) %>% unique())) %>%
    ggplot(aes(x=Description, y=EnrichmentRatio, fill=neglog10p, label=Description)) +
    geom_col() + 
    geom_text(y=max_enrichment*0.05, hjust = 0, size = 3.5)+
    scale_y_continuous(limits = c(0, max_enrichment))+
    theme_bw(base_line_size = 0.75/.px, base_size = 10) +
    scale_fill_continuous(type = "gradient",
                          high = "#4A7AB7", low = "#C2E3EE",
                          space = "Lab", na.value = "grey70", guide = "colourbar",
                          label = ~round(., digits = 1),
                          name = "-Log10(adj p-value)")  +
    coord_flip() +
    ggtitle("ORA: enrichment of downregulated genes") +
    xlab("") +
    ylab("Enrichment Ratio") +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black", linewidth = 0.75/.px),
          axis.text = element_text(color = "black", size = 10),
          legend.text = element_text(color = "black", size = 10),
          plot.title.position = "plot",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          title = element_text(size = 10))
  
  
  # Fixed dimensions based on number of pathways
  plot.pos <-  set_panel_size(plot.pos,
                              # define plot hieght based on the numer of pathways to be plotted
                              height=unit((.ORA.pos %>% filter(qvalue <= cutoff ) %>% nrow())/2, "cm"),
                              width=unit(5, "cm"))
  plot.neg <-  set_panel_size(plot.neg,
                              # define plot hieght based on the numer of pathways to be plotted
                              height=unit((.ORA.neg %>% filter(qvalue <= cutoff ) %>% nrow())/2, "cm"),
                              width=unit(5, "cm"))
  
  
  if (subgroup == "Up") {
    .plot <- grid.arrange(plot.pos, ncol=1, top=grid::textGrob(title),
                          # Add 2 to free up some space for the title
                          heights =c((.ORA.pos %>% filter(qvalue <= cutoff ) %>% nrow()/2) + 3))
  } else if (subgroup == "Down") {
    .plot <- grid.arrange(plot.neg, ncol=1, top=grid::textGrob(title),
                          # Add 2 to free up some space for the title
                          heights =c((.ORA.neg %>% filter(qvalue <= cutoff ) %>% nrow()/2) + 3))
  } else if (subgroup == "all") {
    # plot the plot
    .plot <- grid.arrange(plot.pos, plot.neg, ncol=1, top=grid::textGrob(title),
                          # Add 2 to free up some space for the title
                          heights =c((.ORA.pos %>% filter(qvalue <= cutoff ) %>% nrow()/2) + 3,
                                     .ORA.neg %>% filter(qvalue <= cutoff ) %>% nrow()/2) + 2)
  }
  plot(.plot)
  # return a list, first item is the plot, second are the dimensions for saving it 
  list("plot" = .plot,
       "height" = case_when(
         # Add some margins for axis and others
         subgroup == "pos" ~ (.ORA.pos %>% filter(qvalue <= cutoff) %>% nrow()/2 + 3),
         subgroup == "neg" ~ (.ORA.neg %>% filter(qvalue <= cutoff) %>% nrow()/2 + 3),
         subgroup == "all" ~ (.ORA %>% filter(qvalue <= cutoff) %>% nrow()/2 + 6)
       ))
}


runORA <- function(result, annotation, title = "", cutoff = 0.05, log2FC_threshold = 1, padj_threshold = 0.05, plot = FALSE, toDataFrame = TRUE, custom_annotation=F, ...){
  if(custom_annotation) {
    # custom t2g object to be passed into annotatation argument
    selected_t2g <- annotation
  } else {
    selected_t2g <- # extract t2g of interest based on string passed into annotation
      t2g %>% filter(gs_subcat == annotation) %>% select(gs_name, ensembl_gene)
  }
  
  .em_pos <- enricher(result %>% filter(log2FoldChange > log2FC_threshold, padj < padj_threshold) %>% pull("ensembl_gene_id"),
                      TERM2GENE = selected_t2g,
                      pAdjustMethod = "fdr",
                      minGSSize = 5, maxGSSize = 500,
                      pvalueCutoff = 1,
                      qvalueCutoff = 1
  )
  
  .em_pos.summary <- as.data.frame(.em_pos) %>% mutate(direction = "Up")
  
  .em_neg <- enricher(result %>% filter(log2FoldChange < -log2FC_threshold, padj < padj_threshold) %>% pull("ensembl_gene_id"),
                      TERM2GENE = selected_t2g,
                      minGSSize = 5, maxGSSize = 500,
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1
  )
  
  .em_neg.summary <- as.data.frame(.em_neg)  %>% mutate(direction = "Down")
  
  .em.summary <- rbind(.em_pos.summary,
                       .em_neg.summary) %>%
    rowwise() %>%
    mutate(EnrichmentRatio = eval(parse(text = GeneRatio))/eval(parse(text = BgRatio)))
  
  # If not specified, use the object name as plot title
  title <- ifelse(title == "", deparse(substitute(result)), title)
  # PLot
  if (plot == TRUE) {
    plot(plotORA(.em.summary, title = title, cutoff = cutoff, ...))
  }
  
  .em.summary
}


# Plot heatmap -----------------------
plotHeatmap <- function(vsd, metadata, diff_exp_result, selected_genes,
                        x_axis_var, grouping_var,
                        title = "", hide_not_expressed = TRUE,
                        plotting_value = c("zscore", "centered_vst", "vst")) {
  
  # ** prepare Heatmap input data **
  # preserve only valid gene symbols and throw a warning for faulty ones
  selected_genes_valid <- selected_genes[selected_genes %in% ensembl2symbol$external_gene_name]
  if (length(selected_genes_valid) < length(unique(selected_genes))) {
    warning(paste0("The following genes were not found in the ensembl2symbol mapping: ",
                   paste(base::setdiff(selected_genes, selected_genes_valid), collapse = ", ")))
  }
  
  selected_ensembl_genes <- ensembl2symbol %>%
    filter(external_gene_name %in% selected_genes_valid) %>%
    filter(!is.na(external_gene_name)) %>%
    # preserve only genes that are present in the VST data
    filter(ensembl_gene_id %in% rownames(vsd)) %>%
    pull(ensembl_gene_id)
  
  # Filter VST data for selected genes
  vsd_filtered <- vsd[selected_ensembl_genes, ]
  
  # rank genes from highest to lowest differentially expressed
  gene_order <- diff_exp_result %>%
    dplyr::filter(ensembl_gene_id %in% selected_ensembl_genes) %>%
    # if a heatmap gene is not present, add it to the list assuming log2FC = 1, pvalue = 1 for ranking purpose
    bind_rows(
      data.frame(ensembl_gene_id = base::setdiff(selected_ensembl_genes, diff_exp_result$ensembl_gene_id),
                 log2FoldChange = 1, pvalue = 1) %>%
        left_join(ensembl2symbol)
    ) %>%
    arrange(desc(log2FoldChange * -log10(pvalue))) %>%
    select(external_gene_name, ensembl_gene_id)
  
  # convert to long format
  heatmap_data <- vsd_filtered %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    gather(key = "sample_name", value = "vst", -ensembl_gene_id) %>%
    right_join(ensembl2symbol, .) %>%
    left_join(metadata) %>%
    # compute z-score and meaen centering
    group_by(external_gene_name, ensembl_gene_id) %>%
    mutate(zscore = scale(vst, center = TRUE, scale = TRUE)[,1],
           centered_vst = scale(vst, center = TRUE, scale = FALSE)[,1]) %>%
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = gene_order$ensembl_gene_id)) %>%
    arrange(ensembl_gene_id) %>%
    # annotate genes with 2 ensemble_gene_id matching the same hgnc symbol
    group_by(external_gene_name, sample_label) %>%
    mutate(duplicated = n() > 1,
           gene_label = ifelse(duplicated, 
                               paste0(external_gene_name, " (", str_remove(ensembl_gene_id, ("(?<=ENS(MUS)?G)0+")), ")"), 
                               external_gene_name)) %>%
    ungroup() %>%
    mutate(gene_label = fct_inorder(gene_label))
  
  # remove genes with 0 counts (constant vst) whose z_score is NaN
  if (hide_not_expressed) {
    heatmap_data <- heatmap_data %>% filter(!is.nan(zscore))
  }
  
  # ** prepare plotting aesthetics **
  y_label =  recode(plotting_value,
                    "zscore" = "z-score (VST)",
                    "centered_vst" = "mean centered VST",
                    "vst" = "VST")
  
  # measured to achieve an aspect ratio = 1
  n_rows = length(unique(heatmap_data$external_gene_name))
  n_cols = length(unique(
    interaction(heatmap_data[[x_axis_var]], heatmap_data[[grouping_var]])
  ))
  
  # ** plot heatmap **
  heatmap <- ggplot(heatmap_data, aes(x=.data[[x_axis_var]], y=gene_label)
  ) +
    geom_tile(aes(fill=.data[[plotting_value]])) +
    # pick diverging or sequential palette based on data type
    {  if (plotting_value == "vst") {
      scale_fill_continuous_sequential(palette = "inferno", name=y_label,
                                       rev = F,
                                       labels = label_number(drop0trailing = TRUE),
                                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
    } else {
      scale_fill_distiller(type = "div", palette = "RdBu", name=y_label,
                           labels = label_number(drop0trailing = TRUE),
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
    }
    } +
    scale_x_discrete(expand=c(0, 0)) + 
    scale_y_discrete(expand=c(0, 0)) +
    theme_bw(base_line_size = 0.75/.px, base_size = 10) +
    theme(panel.grid          = element_blank(),
          plot.title.position = "plot",
          axis.title.x        = element_blank(),
          axis.ticks          = element_line(linewidth=0.75/.px, color = "black"),
          panel.border        = element_rect(linewidth=0.75/.px, color = "black"),
          panel.spacing.x     = unit(0, "lines"), # remove panel spacing
          legend.position     = "bottom",
          strip.background    = element_blank(),
          strip.text          = element_text(size = 10, color = "black"),
          strip.clip          = "off",
          axis.text           = element_text(size = 10, color = "black"),
          axis.text.y           = element_text(size = 8, color = "black"),
          axis.text.x         = element_text(angle = 90, hjust = 1, vjust = 0.5),
    ) +
    ylab(y_label)+
    facet_grid(~.data[[grouping_var]], scales= "free", space = "free") +
    ggtitle(title) +
    coord_cartesian(clip = "off") + # remove panel clipping mask
    # set  plot dimensions to obtain cubic tiles
    plot_layout(heights = unit(n_rows/3, "cm"), widths = unit(n_cols/3, "cm")) 
  
  heatmap
}
