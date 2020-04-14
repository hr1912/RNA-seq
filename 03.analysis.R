setwd('/workspace/RNA-seq')

library(magrittr)
library(ggplot2)
library(org.Hs.eg.db)

id2name <- readr::read_delim("group1_r1.genelevel_table.txt", delim = "\t") %>% 
    dplyr::select(`Gene Name`, `Gene ID`) %>% dplyr::rename(Symbol = `Gene Name`,  Id = `Gene ID`)


library('DESeq2')

countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))

colData <- as.matrix(c("group1", "group1", "group1", "group2", "group2","group2"))

colnames(colData) <- "group"


rownames(colData) <- colnames(countData)


dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData_tmp, design = ~ group)

dds <- DESeq(dds)
res <- results(dds)

# cutoff with adjusted p value .05 and absolute fold change 2

res %>% tibble::as_tibble(rownames = "Id") %>% dplyr::left_join(id2name) %>% 
    dplyr::filter(padj < .05 & abs(log2FoldChange) >= log2(2)) -> res_filterd

 res %>% tibble::as_tibble(rownames = "Id") %>% dplyr::left_join(id2name) %>% 
     dplyr::filter(!is.na(padj)) -> res_raw


# volcano plot 
 
 res_raw %>% 
     dplyr::filter(!grepl("_", Id)) %>% 
     #dplyr::mutate(qval = ifelse(qval ==0, 1e-322, qval)) %>% 
     dplyr::mutate(diff_mark_color = ifelse(abs(log2FoldChange) >= 1 & padj <.05, "grey50", "lightgrey")) %>% 
     dplyr::mutate(diff_mark_color = ifelse(Symbol %in% (res_filterd %>% dplyr::arrange(log2FoldChange, padj) %>% 
                                                             dplyr::slice(1:20) %>% dplyr::pull(Symbol)), "blue", diff_mark_color)) %>% 
     dplyr::mutate(diff_mark_color = ifelse(Symbol %in% (res_filterd %>% dplyr::arrange(desc(log2FoldChange), padj) %>% 
                                                             dplyr::slice(1:20) %>% dplyr::pull(Symbol)), "red", diff_mark_color)) ->
     res_raw_tbl


 res_raw_tbl %>% 
     ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
     geom_point(size = .5, color = res_raw_tbl$diff_mark_color) +
     ggrepel::geom_text_repel(data=  res_raw_tbl %>% 
                                  dplyr::filter(Symbol %in% (res_filterd %>% dplyr::arrange(log2FoldChange, padj) %>% 
                                                                 dplyr::slice(1:20) %>% dplyr::pull(Symbol) %>% 
                                                                 c(res_filterd %>% 
                                                                       dplyr::arrange(desc(log2FoldChange), padj) %>% 
                                                                       dplyr::slice(1:20) %>% 
                                                                       dplyr::pull(Symbol)))),
                              aes(label=Symbol), size = 3, force = .5,
                              point.padding = 1) +
     #ggsci::scale_color_nejm() +                           
     xlab("log2(fold change)") + ylab("-log10(adjusted p value)") +
     ggtitle("Volcano Plot") +
     theme_classic(base_size = 12) -> p 


library(clusterProfiler)
library(org.Hs.eg.db)
#library(DOSE)

# code for generate enrichment plot of GO

upGenelst <- res_filterd %>% dplyr::filter(log2FoldChange > 0 ) %>% dplyr::pull(Symbol)

upGo <- clusterProfiler::enrichGO(gene  = upGenelst,
                                  keyType = "SYMBOL",
                                  OrgDb = "org.Hs.eg.db",
                                  ont = "BP",
                                  pAdjustMethod = "bonferroni",
                                  pvalueCutoff = .05,
                                  qvalueCutoff = .1)

library(enrichplot)

upGo@result %>% head(5) %>% 
    dplyr::mutate(ID = factor(ID, levels = rev(ID))) %>% 
    ggplot(aes(x= ID,y = -log10(pvalue))) +
    geom_col(col = "darkred", fill = "white") +
    geom_label(aes(label = Description), size =4, position = position_stack(vjust = 0.5)) +
    coord_flip() +
    ggtitle("Up regulation GO") +
    theme_classic(base_size = 16) -> p1


downGenelst <- res_filterd %>% dplyr::filter(log2FoldChange < 0 ) %>% dplyr::pull(Symbol)

downGo <- clusterProfiler::enrichGO(gene  = downGenelst,
                                    keyType = "SYMBOL",
                                    OrgDb = "org.Hs.eg.db",
                                    ont = "BP",
                                    pAdjustMethod = "bonferroni",
                                    pvalueCutoff = .05,
                                    qvalueCutoff = .1)


downGo@result%>% head(5) %>% 
    dplyr::mutate(ID = factor(ID, levels = rev(ID))) %>% 
    ggplot(aes(x= ID,y = -log10(pvalue))) +
    geom_col(col = "darkblue", fill = "white") +
    geom_label(aes(label = Description), size =4, position = position_stack(vjust = 0.5)) +
    coord_flip() +
    ggtitle("Down regulation GO") +
    theme_classic(base_size = 16) -> p2


# code for generating heatmap

colData %>% 
    scale() %>% 
    tibble::as.tibble(rownames = "Id") %>% 
    dplyr::left_join(id2name) %>% 
    dplyr::filter(Symbol %in% res_filterd$Symbol)-> countData_df

res_filterd %>% dplyr::arrange(desc(log2FoldChange)) %>% 
    dplyr::pull(Symbol) %>% unique() -> genes_sorted

countData_df %>% dplyr::select(-Symbol, -Id) %>% 
    as.matrix() -> countData_mat

rownames(countData_mat) <- countData_df$Symbol

countData_mat_scaled <- t(countData_mat) %>% scale() %>% t()

countData_mat_scaled

breakList <- seq(-3, 3, by = .2)


pheatmap::pheatmap(
	countData_mat_scaled,
    countData_mat_scaled[genes_sorted,], 
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(length(breakList)),
    breaks = breakList, 
    cluster_rows = F,
    clustering_method = "complete",
    #                   clustering_method = "centroid", 
    show_rownames = F, fontsize = 16) -> p3

# code for plot several figures together

library(gridExtra)

pdf("Heatmap_with_GO.pdf", width = 9, height = 7)

grid.arrange(p3[[4]], p1, p2, ncol = 2,
             layout_matrix = cbind(c(1,1), c(2,3)),
             widths = c(1.5, 4))

dev.off()

