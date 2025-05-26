library(msigdbr)
library(clusterProfiler)
library(patchwork)
library(ggrepel)
library(egg)
library(gridExtra)
library(tidyverse)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")

###############################
#          OPTIONS            #
###############################

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
VolcanoPlot <- function(result, title = element_blank(), thrLog2FC = 1, thrPadj = 0.05, pctLabel = 0.005, protein_coding_label_only = FALSE) {
  
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
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 1),
          legend.position = "none") +
    xlim(x_limits) +
    scale_color_manual(values = c("dodgerblue", "brown1"), na.value = "grey") +
    xlab(expression(log["2"](Fold~Change))) + ylab(expression(-log["10"](Adjusted~p~Value))) + 
    ggtitle(title)
}

## Prepare imput data for GSEA
PrepareForGSEA <- function(result, .na.rm = FALSE) {
  
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
    theme_bw() +
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

plotGSEAwide <- function(.GSEA, title = "", cutoff = 0.05, subgroup = "all") {
  
  if (subgroup == "pos") {
    .GSEA <- .GSEA %>% filter(NES > 0)
  } else if (subgroup == "neg") {
    .GSEA <- .GSEA %>% filter(NES < 0)
  }
  
  .GSEA %>%
    filter(qvalues <= cutoff) %>%
    mutate(Description = factor(Description, levels = 
                                  .GSEA %>% arrange(desc(NES)) %>% pull(Description))) %>%
    ggplot(aes(x=Description, y=NES, fill=qvalues)) +
    geom_col() + 
    theme_bw() +
    scale_fill_continuous(type = "gradient",
                          low = "#E74C3C", high = "#F1C40F",
                          space = "Lab", na.value = "grey70", guide = "colourbar")  +
    #   coord_flip() +
    ggtitle(paste0("GSEA: ", title)) +
    xlab("") +
    ylab("Normalized Enrichment Score") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin=unit(c(0.2,0.2,0.2,2),"cm"),
          plot.title.position = "plot")
}



## Run GSEA
runGSEA <- function(result, annotation, title = "", cutoff = 0.05, plot = FALSE, toDataFrame = TRUE, custom_annotation=F, ...) {
  if(custom_annotation) {
    # custom t2g object to be passed into annotatation argument
    selected_t2g <- annotation
  } else {
    selected_t2g <- # extract t2g of interest based on string passed into annotation
      t2g %>% filter(gs_subcat == annotation) %>% select(gs_name, ensembl_gene)
  }
    
  .em <- GSEA(PrepareForGSEA(result),
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

# Define custom GSEA plot function ------------
custom_gseaplot <- function(res, pathway) {
  statistics <- res %>% as.data.frame() %>% filter(Description == pathway)
  gseaplot(res, geneSetID = pathway, by = "runningScore", title = pathway, color.line = "black") + 
    theme_minimal() + 
    ylab("Enrichment Score") +
    xlab("")+
    labs(caption = paste0("NES: ", round(statistics$NES,2),", p-value: ", round(statistics$pvalue,4), ", q-value: ",  round(statistics$qvalues,4)))
}

# Plot heatmap -----------------------
plotHEATMAP <- function(zscore, selected_genes, title, ylab) {
  heatmap.input <- zscore %>%
    filter(external_gene_name %in% selected_genes,
           !is.na(zscore))
  
  heatmap.input <- heatmap.input %>%
    mutate(external_gene_name = factor(external_gene_name,
                                       levels=result %>% arrange(log2FoldChange*-log10(pvalue)) %>% filter(external_gene_name %in% selected_genes) %>% pull("external_gene_name") %>% unique()))
  
  heatmap <- heatmap.input %>%
    # remove cell line name
  #  mutate(sample =  str_replace_all(sample, "_", " ") %>% str_remove("^[^\\s]+")) %>%
    ggplot(aes(x=sample, y=external_gene_name, fill=zscore)) +
    geom_tile()+
    scale_x_discrete(expand=c(0, 0)) + 
    scale_y_discrete(expand=c(0, 0))+
    scale_fill_distiller(type = "div",palette = 5, name="z-score(log2(CPM+1))",
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    theme(legend.position = "none")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot",
          panel.border = element_rect(linewidth=0.75),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    ylab(ylab)+
    xlab("") +
    ggtitle(title)
  
  heatmap <-  set_panel_size(heatmap,
                             # define plot hieght based on the numer of pathways to be plotted
                             height=unit((selected_genes %in% heatmap.input$external_gene_name %>% sum())/2, "cm"),
                             width=unit(8, "cm"))
  
  heatmap
}

## Plot GSEA table
plotGSEAcompact <- function(.GSEA, title = "", cutoff = 0.05, subgroup = "all", truncate_label_at = 45) {
  .GSEA <- .GSEA %>% 
    filter(qvalues <= cutoff) %>%
    # Truncate long strings
    mutate(Description = str_trunc(Description, width = truncate_label_at, side = "right", ellipsis = ".."))
  
  # use it to set a common range scale
  max_NES <- .GSEA %>% filter(qvalues <= cutoff) %>% top_n(1, abs(NES)) %>% pull("NES") %>% abs() %>% unique()
  
  .GSEA.pos = .GSEA %>% filter(NES > 0)
  .GSEA.neg = .GSEA %>% filter(NES < 0)
  
  plot.pos <- .GSEA.pos %>%
    mutate(neglog10p = -log10(qvalues),
           # set in gray non significant pathway
           neglog10p = ifelse(neglog10p < -log10(0.05), as.double(NA), neglog10p),
           Description = factor(Description, levels = 
                                  .GSEA.pos %>% arrange(NES) %>% pull(Description))) %>%
    ggplot(aes(x=Description, y=NES, fill=neglog10p, label=Description)) +
    geom_col() + 
    geom_text(y=max_NES*0.05, hjust = 0, size = 3)+
    scale_y_continuous(limits = c(0, max_NES))+
    theme_bw(8) +
    scale_fill_continuous(type = "gradient",
                          high = "#DC3D2D", low = "#FED98B",
                          space = "Lab", na.value = "grey70",
                          label =  ~round(., digits = 1),
                          guide = "colourbar", name = "-Log10(adj p-value)")  +
    coord_flip() +
    ggtitle("GSEA: positive enrichment") +
    xlab("") +
    ylab("Normalized Enrichment Score") +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.5, color="black", fill=NA),
          plot.title.position = "plot",
          axis.text = element_text(color = "black", size = 8),
          axis.text.y = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.ticks.y = element_blank(),
          title = element_text(size = 9))
  
  plot.neg <- .GSEA.neg %>%
    mutate(neglog10p = -log10(qvalues),
           neglog10p = ifelse(neglog10p < -log10(0.05), as.double(NA), neglog10p),
           Description = factor(Description, levels = 
                                  .GSEA.neg %>% arrange(desc(NES)) %>% pull(Description))) %>%
    ggplot(aes(x=Description, y=NES, fill=neglog10p, label=Description)) +
    geom_col() + 
    geom_text(y=max_NES*-0.05, hjust = 1, size = 3)+
    scale_y_continuous(limits = c(-max_NES, 0))+
    theme_bw(8) +
    scale_fill_continuous(type = "gradient",
                          high = "#4A7AB7", low = "#C2E3EE",
                          space = "Lab", na.value = "grey70", guide = "colourbar",
                          label = ~round(., digits = 1),
                          name = "-Log10(adj p-value)")  +
    coord_flip() +
    ggtitle("GSEA: negative enrichment") +
    xlab("") +
    ylab("Normalized Enrichment Score") +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.5, color="black", fill=NA),
          plot.title.position = "plot",
          axis.text = element_text(color = "black", size = 8),
          axis.text.y = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.ticks.y = element_blank(),
          title = element_text(size = 9)
          )
  
  # Fixed dimensions based on number of pathways
    plot.pos <-  set_panel_size(plot.pos,
                            # define plot hieght based on the numer of pathways to be plotted
                            height=unit((.GSEA.pos %>% filter(qvalues <= cutoff ) %>% nrow())/1.5, "cm"),
                            width=unit(12, "cm"))
    plot.neg <-  set_panel_size(plot.neg,
                                # define plot hieght based on the numer of pathways to be plotted
                                height=unit((.GSEA.neg %>% filter(qvalues <= cutoff ) %>% nrow())/1.5, "cm"),
                                width=unit(12, "cm"))
    
    
    if (subgroup == "pos") {
      .plot <- grid.arrange(plot.pos, ncol=1, top=grid::textGrob(title),
                            # Add 2 to free up some space for the title
                            heights =c((.GSEA.pos %>% filter(qvalues <= cutoff ) %>% nrow()/1.5) + 3))
    } else if (subgroup == "neg") {
      .plot <- grid.arrange(plot.neg, ncol=1, top=grid::textGrob(title),
                            # Add 2 to free up some space for the title
                            heights =c((.GSEA.neg %>% filter(qvalues <= cutoff ) %>% nrow()/1.5) + 3))
    } else if (subgroup == "all") {
      # plot the plot
      .plot <- grid.arrange(plot.pos, plot.neg, ncol=1, top=grid::textGrob(title),
                            # Add 2 to free up some space for the title
                            heights =c((.GSEA.pos %>% filter(qvalues <= cutoff ) %>% nrow()/1.5) + 3,
                                       .GSEA.neg %>% filter(qvalues <= cutoff ) %>% nrow()/1.5) + 2)
    }
    plot(.plot)
    # return a list, first item is the plot, second are the dimensions for saving it 
    list("plot" = .plot,
         "height" = case_when(
           # Add some margins for axis and others
           subgroup == "pos" ~ (.GSEA.pos %>% filter(qvalues <= cutoff) %>% nrow()/1.5 + 3),
           subgroup == "neg" ~ (.GSEA.neg %>% filter(qvalues <= cutoff) %>% nrow()/1.5 + 3),
           subgroup == "all" ~ (.GSEA %>% filter(qvalues <= cutoff) %>% nrow()/1.5 + 6)
         ))
}

# ORA -----------------------------------

## Plot ORA table
plotORA <- function(.ORA, title = "", cutoff = 0.05, subgroup = "all", truncate_label_at = 45) {
  .ORA <- .ORA %>% 
    filter(qvalue <= cutoff) %>%
    # Truncate long strings
    mutate(Description = str_trunc(Description, width = truncate_label_at, side = "right", ellipsis = ".."))
  
  # use it to set a common range scale
  max_ratio <- .ORA %>% filter(qvalue <= cutoff) %>% pull("ratio") %>% max()
  
  .ORA.pos = .ORA %>% filter(direction == "pos")
  .ORA.neg = .ORA %>% filter(direction == "neg")
  
  plot.pos <- .ORA.pos %>%
    mutate(neglog10p = -log10(qvalue),
           Description = factor(Description, levels = 
                                  .ORA.pos %>% arrange(ratio) %>% pull(Description) %>% unique())) %>%
    ggplot(aes(x=Description, y=ratio, fill=neglog10p, label=Description)) +
    geom_col() + 
    geom_text(y=max_ratio*0.05, hjust = 0, size = 3.5)+
    scale_y_continuous(limits = c(0, max_ratio), labels = scales::label_percent())+
    theme_bw() +
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
          plot.title.position = "plot",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          title = element_text(size = 9))
  
  plot.neg <- .ORA.neg %>%
    mutate(neglog10p = -log10(qvalue),
           Description = factor(Description, levels = 
                                  .ORA.neg %>% arrange(ratio) %>% pull(Description) %>% unique())) %>%
    ggplot(aes(x=Description, y=ratio, fill=neglog10p, label=Description)) +
    geom_col() + 
    geom_text(y=max_ratio*0.05, hjust = 0, size = 3.5)+
    scale_y_continuous(limits = c(0, max_ratio), labels = scales::label_percent())+
    theme_bw() +
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
          plot.title.position = "plot",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          title = element_text(size = 9))
  
  
  # Fixed dimensions based on number of pathways
  plot.pos <-  set_panel_size(plot.pos,
                              # define plot hieght based on the numer of pathways to be plotted
                              height=unit((.ORA.pos %>% filter(qvalue <= cutoff ) %>% nrow())/2, "cm"),
                              width=unit(5, "cm"))
  plot.neg <-  set_panel_size(plot.neg,
                              # define plot hieght based on the numer of pathways to be plotted
                              height=unit((.ORA.neg %>% filter(qvalue <= cutoff ) %>% nrow())/2, "cm"),
                              width=unit(5, "cm"))
  
  
  if (subgroup == "pos") {
    .plot <- grid.arrange(plot.pos, ncol=1, top=grid::textGrob(title),
                          # Add 2 to free up some space for the title
                          heights =c((.ORA.pos %>% filter(qvalue <= cutoff ) %>% nrow()/2) + 3))
  } else if (subgroup == "neg") {
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
  
  .em_pos.summary <- as.data.frame(.em_pos) %>% mutate(direction = "pos")
  
  .em_neg <- enricher(result %>% filter(log2FoldChange < -log2FC_threshold, padj < padj_threshold) %>% pull("ensembl_gene_id"),
                      TERM2GENE = selected_t2g,
                      minGSSize = 5, maxGSSize = 500,
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1
  )
  
  .em_neg.summary <- as.data.frame(.em_neg)  %>% mutate(direction = "neg")
  
  .em.summary <- rbind(.em_pos.summary,
                       .em_neg.summary) %>%
    mutate(ratio = as.numeric(str_extract(GeneRatio, "^[0-9]+"))/as.numeric(str_extract(BgRatio, "^[0-9]+")))
  
  # If not specified, use the object name as plot title
  title <- ifelse(title == "", deparse(substitute(result)), title)
  # PLot
  if (plot == TRUE) {
    plot(plotORA(.em.summary, title = title, cutoff = cutoff, ...))
  }
  
  .em.summary
}

