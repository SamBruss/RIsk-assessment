# Install required packages
library(RCy3)
library(dplyr)
library(rWikiPathways)
library(readxl)
## In your current directory, create a file name "File_Data" to put each data set
## add the linkset file for extension in a file name "Linkset"
  
cytoscapePing()
installApp('WikiPathways')

## In this string store the name of each data set

omics_data <- c("GSE144775.PFOA.0.02.D1.top.table.tsv", 
                "GSE144775.PFOA.0.1.D1.top.table.tsv", 
                "GSE144775.PFOA.0.2.D1.top.table.tsv",
                "GSE144775.PFOA.1.D1.top.table.tsv",
                "GSE144775.PFOA.10.D1.top.table.tsv",
                "GSE144775.PFOA.100.D1.top.table.tsv")

# create excel sheet to display Key Enrichment score and set the name of the finale file
KE_SCORE <- createWorkbook()
xls_name <- "PFAS_RISK"
WP = "WP5464"
V = 1
# create for loop for every data sheet 
for (data_files in omics_data) {

# now import AOP as a network using the file directory 
  commandsRun(paste0('wikipathways import-as-network id=', WP))

### Before extension a new columns is created with the WPID only for extension mapping 
# separate WP ID column
  nodetable <- getTableColumns("node",c("SUID","XrefId"))
# Initialize WPID column to store results
  nodetable$WPID <- NA
#filter and transfer the data to WPID
  nodetable <- nodetable %>%
  mutate(WPID = ifelse(grepl("[a-zA-Z]", XrefId), XrefId, WPID),
  XrefId = ifelse(grepl("[a-zA-Z]", XrefId), NA, XrefId))
  print(nodetable)

# add a key WP ID table to cytoscape
  loadTableData(nodetable, data.key.column = "SUID", table = "node", table.key.column = "SUID")

# extend pathway using WPID as attribute
  linkset_dir <- paste0(getwd(),"/Linkset/WikiPathways_Hsa_20230910.xgmml")
  CTLextend.cmd = paste('cytargetlinker extend idAttribute="WPID" linkSetFiles="', linkset_dir, '" network=current direction=TARGETS')
  commandsRun(CTLextend.cmd)
  layoutNetwork()

# Add omics data and set the name of the network
  data_name <- paste(data_files)
  genedata_dir <- paste0(getwd(), "/Final_Data/", data_name)
  netwok_name <- paste(data_name)
  gene_data <- read.delim(genedata_dir)
  loadTableData(gene_data, data.key.column = "GeneID", table = "node", table.key.column = "CTL.GeneID")
  renameNetwork(data_name)

# before adding style, a "Type" column is create for better mapping 
  mapping.table <- getTableColumns("node",c("SUID","Type"))
  mapping.table$Type[is.na(mapping.table$Type)] <- 'Gene'
  mapping.table$Type[mapping.table$Type == 'Unknown'] <- 'Event'
  # the column is added to cytoscape
  loadTableData(mapping.table, data.key.column = "SUID", table = "node", table.key.column = "SUID")

# add style 
createVisualStyle("Gene_vis")
setVisualStyle("Gene_vis")

# Label
setNodeLabelMapping(
  table.column = 'name',
  style.name = 'Gene_vis'
)

# Shape for different type of nodes 

setNodeShapeMapping(
  table.column = 'Type',
  table.column.values = c("Event", "Gene", "Pathway"),
  shapes = c("RECTANGLE", "ELLIPSE", "DIAMOND"),
  style.name = "Gene_vis"
)

# set color mapping

setNodeColorMapping(
  table.column = 'log2FoldChange',
  table.column.values = c('-3.0', '0.0', '3.0'),
  colors = c('#3E8CE0', '#FFFFFF', '#FE6565'),
  mapping.type = "c",
  style.name = "Gene_vis",
  default.color = "#FFFFFF"
)

# set color mapping for key events and pathway 
  setNodeBorderColorMapping(
    table.column = 'Type',
    table.column.values = c("Event", "Pathway"),
    colors = c('#000000', '#1CBF23'),
    mapping.type = "d",
    style.name = "Gene_vis"
)

# filter gene with no values

gene_filter <- getTableColumns(table = "node",
                columns = c("SUID", "pvalue","Type"))
gene_filter <- gene_filter %>%
  filter(!(Type %in% c("Event", "Pathway")))

# check if the gene have pvalue
gene_filter <- gene_filter %>%
  mutate(is_gene_pvalue = grepl("gene", Type, ignore.case = TRUE) & is.na(pvalue))
# check for significant genes
gene_filter <- gene_filter %>%
  mutate(significant = ifelse(pvalue < 0.05, TRUE, NA))
gene_filter <- gene_filter %>%
  select(-pvalue)
# Load mapping table in cytoscape for style
loadTableData(data = gene_filter, data.key.column = "SUID",
              table.key.column = "SUID",
              table = "node")

# Set style for significant genes
setNodeBorderWidthMapping(table.column = "significant",
                          table.column.values = "true",
                          widths = "5",
                          mapping.type = "d",
                          style.name = "Gene_vis")
createColumnFilter(filter.name = "no_values_genes",
                   column = "is_gene_pvalue",
                   criterion = TRUE,
                   predicate = "IS")
deleteSelectedNodes()
clearSelection()


############################
############################
# KE Enrichment calculation
# Each gene linker in the node table (pathway, GO annotation, KEGG, etc) should be name "Pathway" in the "Type" column

value_file <- data.frame(b = NA,
                         B = NA,
                         n = NA,
                         N = NA,
                         row.names = "values")

# add the population number value (N)
createColumnFilter(filter.name = "allgenes", 
                   column = "Type",
                   criterion = "Gene",
                   predicate = "IS")

value_file$N <- getSelectedNodeCount()

clearSelection()

# add the population success number value (n)
createColumnFilter(filter.name = "siggenes",
                   column = "pvalue",
                   criterion = 0.05,
                   predicate = "LESS_THAN")

value_file$n <- getSelectedNodeCount()

clearSelection()


#now lets find the b and B value for each separate pathway

#get the table with name, pvalue, SUID and Type to later filter the gene
pvalue_file <- getTableColumns("node", c("SUID", "pvalue", "name", "Type"))

# retrieve the pathway or Gene Ontology nodes

createColumnFilter(filter.name = "Pathway",
                   column = "Type",
                   criterion = "Pathway",
                   predicate = "IS")
pathway_analysis <- data.frame(Pathway = getSelectedNodes())

clearSelection()

# a for loop that retreive the b and B values of each pathway
# first the variable X is set to one which will define the row to store KE enrichment score of the corresponding pathway 

  for (y in pathway_analysis$Pathway) {
  
  selectNodes(by.col = "name", sprintf('%s', y))
  connected_gene <- data.frame(name = getFirstNeighbors())
  
  pvalue_connected_gene <- connected_gene %>%
    left_join(pvalue_file, by = "name")
  pvalue_connected_gene <- pvalue_connected_gene %>% filter(!is.na(pvalue))
  
  pvalue_connected_gene <- pvalue_connected_gene %>%
    mutate(significant = pvalue < 0.05)
  
  # add sample size number value (B)
  value_file$B <- nrow(pvalue_connected_gene)
  
  # add the number of success in the sample (b)
  if (sum(pvalue_connected_gene$significant == TRUE) == 0) {
    value_file$b <- sum(pvalue_connected_gene$significant == TRUE) + 1
  } else {
    value_file$b <- sum(pvalue_connected_gene$significant == TRUE)
  }
  
  # calculate and add KE enrichment score to the pathway table 
  value_file$KE_ENR <- (value_file$b/value_file$B)/(value_file$n/value_file$N)
  print(paste("KE enrichment score of", y, "=", value_file$KE_ENR, "with b = ", value_file$b,
              "B = ", value_file$B, "n = ", value_file$n, "and N = ", value_file$N))
  pathway_analysis$KE_ENR[X] <- print(value_file$KE_ENR)
  
  # calculate Hypergeometric pvalue
  pathway_analysis$HyperG[X] <- dhyper(x = value_file$b,
                                       m = value_file$n ,
                                       n = value_file$N - value_file$n,
                                       k = value_file$B)
  
  pathway_analysis <- pathway_analysis %>%
    mutate(KE_ENR = ifelse(`HyperG` < 0.05, paste0(KE_ENR, "*"), KE_ENR))
  
  
  X=X+1
  clearSelection()
  

  }

addWorksheet(KE_SCORE, sheetName = sprintf("Network%s", V ))
writeData(KE_SCORE, sheet = sprintf("Network%s", V ), x = pathway_analysis)
sheet_name[X] <- (sprintf("Network%s", V ))
V = V+1

layoutNetwork()


}

saveWorkbook(KE_SCORE, sprintf("%s.xls", xls_name), overwrite = TRUE)

openXL(file = sprintf("%s", xls_name))
