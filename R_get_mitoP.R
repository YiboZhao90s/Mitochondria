library(biomaRt)
setwd("C:/Users/yz735/OneDrive - University of Leicester/Fellowship/UKDRI_AI")

mitoGolist <- read.delim("mitoGOlist.txt", stringsAsFactors = F, header = F)
go_id <- mitoGolist[,1]


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
go_data <- getBM(attributes = c("go_id", "name_1006", "ensembl_gene_id", "external_gene_name"),
                 filters = "go",
                 values = go_id,
                 mart = ensembl)
head(go_data)
length(unique(go_data$external_gene_name))



BiocManager::install("ontoProc")
library(ontoProc)
library(GO.db)
ontology_lookup <- select(GO.db, keys = go_id, columns = "ONTOLOGY", keytype = "GOID")

# Split by ontology
bp_terms <- ontology_lookup$GOID[ontology_lookup$ONTOLOGY == "BP"]
cc_terms <- ontology_lookup$GOID[ontology_lookup$ONTOLOGY == "CC"]
mf_terms <- ontology_lookup$GOID[ontology_lookup$ONTOLOGY == "MF"]

paste0(bp_terms, collapse = ",")
paste0(cc_terms, collapse = ",")
paste0(mf_terms, collapse = ",")

write.table(unique(go_data$external_gene_name), "mitoGenelist.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(go_data, "mitoGeneAnnot.txt", row.names = F, quote = F, sep = "\t")

# disease enrichment analysis
#BiocManager::install("DOSE")
library(DOSE)
library(org.Hs.eg.db)

# Example gene symbols
gene_symbols <- unique(go_data$external_gene_name)

# Convert to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
data(geneList)
x <- enrichDGN(gene          = entrez_ids,
               pvalueCutoff  = 0.05,
               pAdjustMethod = "BH",
               universe      = entrez_ids,
               minGSSize     = 5,
               maxGSSize     = 500,
               qvalueCutoff  = 0.05,
               readable      = TRUE)
head(x, 10L)
#BiocManager::install("enrichplot")
library(enrichplot)
dotplot(x, showCategory = 10, title = "Disease Enrichment Dotplot (Top 10)")

write.table(x, "DOSE_res.txt", sep = "\t", quote = F, row.names = F)
