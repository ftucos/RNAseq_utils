library(biomaRt)
library(data.table)
library(tidyverse)
library(glue)

##############################
#           OPTIONS          #
##############################

useCache <- TRUE
cacheDir <- "/Users/tucos/.bioinfoCache"

##############################
#             RUN            #
##############################

# Look for the latest cache files available based on date of data retrieval
biomartCache.human <- list.files(cacheDir, "ensembl2symbol-hsapiens_gene_ensembl.tsv", full.names = T) %>% sort(decreasing = T) %>% head(n=1)
biomartCache.mouse <- list.files(cacheDir, "ensembl2symbol-mmusculus_gene_ensembl.tsv", full.names = T) %>% sort(decreasing = T) %>% head(n=1)
t2gCache.human <- list.files(cacheDir, "t2g-human.tsv.gz", full.names = T) %>% sort(decreasing = T) %>% head(n=1)
t2gCache.mouse <- list.files(cacheDir, "t2g-mouse.tsv.gz", full.names = T) %>% sort(decreasing = T) %>% head(n=1)

date <- Sys.Date()

# specify preferred packages to avoid names conflicts
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")

if (species %in% c("human", "homo sapiens", "hsapiens", "hsa")) {
  
  # if cache is available for human biomart annotation used for ensembl_id to hgnc_symbol conversion, load it 
  if (useCache & !is.na(biomartCache.human)){
    biomartCache <- biomartCache.human
  } else {
    useCache == FALSE
  }
  # if cache is available for human term 2 gene MSigDB annotation for pathway analysis, load it 
  if (useCache & !is.na(t2gCache.human)){
    t2gCache <- t2gCache.human
  }  else {
    useCache == FALSE
  }
  
  # 
  biomart_species = "hsapiens_gene_ensembl"
  kegg_species = "hsa"
  
} else if (species %in% c("mouse", "mus musculus", "mmusculus", "mmu")) {
  if (useCache & !is.na(biomartCache.mouse)){
    biomartCache <- biomartCache.mouse
  } else {
    useCache == FALSE
  }
  if (useCache & !is.na(t2gCache.mouse)){
    t2gCache <- t2gCache.mouse
  }  else {
    useCache == FALSE
  }
  
  biomart_species = "mmusculus_gene_ensembl"
  kegg_species = "mmu"
} else {
  stop(glue("Error: {species} not allowed. Species must be either 'human' or 'mouse'"))
  
}

# Load biomart annotation used for Ensembl Id (ensembl_gene_id) to HGNC Symbol (external_gene_name) conversion -----------------------
if(useCache) {
  # Load cached annotation
  ensembl2symbol = fread(biomartCache)
} else {
  # Load the latest biomart annotation
  mart <- biomaRt::useMart("ensembl", dataset = biomart_species)
  ensembl2symbol <- biomaRt::getBM(mart = mart,
                                   attributes = c("external_gene_name", "ensembl_gene_id",
                                                  "gene_biotype",
                                                  "description")) %>%
    ## Remove source info from the description field
    mutate(description = str_remove(description, "\\s?\\[.*\\]\\s?"))
  # write the file to cache
  write_tsv(ensembl2symbol, paste0(cacheDir, "/", date, "-ensembl2symbol", "-", biomart_species, ".tsv"))
}

# GSEA annotation data -----------------------------------------------------------------
# Load latest KEGG pathway annotation data -----------------------------------
if(useCache) {
  # Load latest KEGG annotation cached
  t2g = fread(t2gCache)
} else {
  # Load KEGG association between pathway ID and gene ID
  pathway2gene = read_tsv(paste0("http://rest.kegg.jp/link/",kegg_species,"/pathway"), col_names = c("pathway_id", "gene_id"))
  # Load KEGG association between pathway ID and pahtway names
  pathway2name = read_tsv(paste0("http://rest.kegg.jp/list/pathway/", kegg_species), col_names = c("pathway_id", "pathway"))
  # Load KEGG association between gene ID and gene names
  gene2name = read_tsv(paste0("http://rest.kegg.jp/list/", kegg_species), col_names = c("gene_id", "biotype", "X3", "gene")) 
  
  gene2name.1 <- gene2name %>%
    # remove genes with no semicolons, they have no official HGNC name (eg. annotated cDNA sequences)
    filter(str_detect(gene, ";")) %>%
    # extract official HGNC name
    mutate(gene = str_extract(gene, "^[^;]+") %>% strsplit(",\\s")) %>%
    # you can have multiple gene names mapping the same id because the mapping is based around orthologues
    unnest(cols = "gene")
  
  # no empty fields
  #table(is.na(gene2name))
  
  # build kegg term to gene table
  kegg_t2g <- pathway2name %>%
    left_join(pathway2gene %>% mutate(pathway_id = str_remove(pathway_id, "path:"))) %>%
    left_join(gene2name.1)
  
  # count number of empty fields
  n_rows_with_missing_genes = sum(is.na(kegg_t2g))
  
  print(glue("Removing {n_rows_with_missing_genes} rows with missing gene names"))
  kegg_t2g <- kegg_t2g %>%
    filter(!is.na(gene))
  
  #table(is.na(kegg_t2g))
  
  # Check how many genes symbols are present in the ensembl2symbol table
  table(kegg_t2g$gene %in% ensembl2symbol$external_gene_name)
  
  kegg_t2g <- kegg_t2g %>%
    mutate(pathway = str_remove(pathway, " - Homo sapiens \\(human\\)") %>% str_remove(" - Mus musculus \\(house mouse\\)"))%>%
    # adapt to the msigdb column names
    select(external_gene_name = gene, gs_name = pathway) %>%
    inner_join(ensembl2symbol %>% select(external_gene_name, ensembl_gene = ensembl_gene_id), by="external_gene_name") %>%
    distinct() %>%
    mutate(gs_subcat = "KEGG_latest") %>%
    dplyr::rename("gene_symbol" = "external_gene_name")
  
  ## Load annotation sets for GSEA
  msigdb <- msigdbr(species = species)
  # clean version of the annotation database
  t2g <- msigdb %>%
    mutate(gs_subcat = ifelse(gs_cat == "H", "HALLMARK", str_remove(gs_subcat, "^CP:"))) %>%
    # select signature of interest
    filter(gs_subcat %in% c("HALLMARK", "KEGG", "REACTOME", "PID", "BP", "WIKIPATHWAYS", "GO:BP")) %>%
    mutate(gs_name = case_when(
      gs_subcat == "HALLMARK" ~ str_remove(gs_name, "^HALLMARK_") %>% str_replace_all("_", " "),
      gs_subcat == "KEGG" ~ gs_description,
      gs_subcat == "REACTOME" ~ gs_description,
      gs_subcat == "PID" ~ gs_description,
      gs_subcat == "WIKIPATHWAYS" ~ gs_description,
      gs_subcat == "GO:BP" ~ str_remove(gs_name, "^GOBP_") %>% str_replace_all("_", " "),
      TRUE ~ gs_name)) %>%
    # add latest kegg annotations
    bind_rows(kegg_t2g)
  
  write_tsv(t2g, paste0(cacheDir, "/", date, "-t2g", "-", species, ".tsv"))
  # compress the file
  system(paste0("gzip ", cacheDir, "/", date, "-t2g", "-", species, ".tsv"))
  
  # remove intermediate files 
  rm(gene2name, kegg_t2g, pathway2gene, pathway2name, biomart_species, kegg_species, date)
}

# Pathway annotation table --------------------------------------------
# Keep only essential annotations of interest

pathway_annotation_table <- full_join(
  t2g %>% filter(gs_subcat == "HALLMARK") %>%
    select(gs_name, ensembl_gene) %>%
    group_by(ensembl_gene) %>%
    summarize(pathways = paste0(gs_name, collapse = "|")) %>%
    dplyr::rename("ensembl_gene_id" = "ensembl_gene",
                  "HALLMARK" = "pathways"),
  t2g %>% filter(gs_subcat == "KEGG_latest") %>%
    group_by(ensembl_gene) %>%
    summarize(pathways = paste0(gs_name, collapse = "|")) %>%
    dplyr::rename("ensembl_gene_id" = "ensembl_gene",
                  "KEGG_latest" = "pathways")) %>%
  full_join(
    t2g %>% filter(gs_subcat == "REACTOME") %>%
      group_by(ensembl_gene) %>%
      summarize(pathways = paste0(gs_name, collapse = "|")) %>%
      dplyr::rename("ensembl_gene_id" = "ensembl_gene",
                    "REACTOME" = "pathways")) %>%
  full_join(
    t2g %>% filter(gs_subcat == "PID") %>%
      group_by(ensembl_gene) %>%
      summarize(pathways = paste0(gs_name, collapse = "|")) %>%
      dplyr::rename("ensembl_gene_id" = "ensembl_gene",
                    "PID" = "pathways"))

