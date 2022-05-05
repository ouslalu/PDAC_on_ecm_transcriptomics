library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(stats)
library(VennDiagram)
library(eulerr)
library(viridis)
library(lubridate)
library(forcats)

read_counts <- read.csv("/Users/myname/Desktop/miapaca on ECM Project (Genomics)/Miapaca_genomics_after_combining_technical_replicates/concatenated_matrix_for_downstream_analysis.csv")
read_counts <- data.frame(read_counts, row.names = "X0")
ecms_of_interest <- c("Fibronectin_1","Fibronectin_2","FBS_1","FBS_2",
                      "B27_1","B27_2","Vitronectin_1","Vitronectin_2", "Vitronectin_3")
read_counts <- read_counts[,ecms_of_interest]

samples<-colnames(read_counts)
group <- c("fibronectin", "fibronectin", "fbs", "fbs", "b27","b27", "vitronectin", 
           "vitronectin","vitronectin")
metadata<-data.frame(sample =samples, group=group, stringsAsFactors = TRUE)
row.names(metadata) <- samples


dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = metadata,
                              design = as.formula("~ group"))


keep <- rowSums(DESeq2::counts(dds))>10
dds <- dds[keep,]

dds <- DESeq(dds)

">>>>>>>>>>>>>>>>>>>>>>>>FBS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"


res_fbs <- results(dds, contrast = c("group", "fbs", "b27"))
res_fbs <- na.omit(res_fbs)
res_fbs_sig <- res_fbs[res_fbs$padj<0.05,]
res_fbs_sig <- res_fbs_sig[abs(res_fbs_sig$log2FoldChange)>0.58,]
upregulated_fbs_sig <- res_fbs_sig[res_fbs_sig$log2FoldChange>0,] 
downregulated_fbs_sig <- res_fbs_sig[res_fbs_sig$log2FoldChange<0,] 
#write.csv(upregulated_fbs_sig, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/fbs_upregulated.csv", quote=FALSE)
#write.csv(downregulated_fbs_sig, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/fbs_downregulated.csv", quote=FALSE)


">>>>>>>>>>>>>>>>>>>>>>>>>Fibronectin>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

res_fibronectin <- results(dds, contrast = c("group", "fibronectin", "b27"))
res_fibronectin <- na.omit(res_fibronectin)
res_fibronectin_sig <- res_fibronectin[res_fibronectin$padj<0.05,]
res_fibronectin_sig <- res_fibronectin_sig[abs(res_fibronectin_sig$log2FoldChange)>0.58,]
upregulated_fibronectin_sig <- res_fibronectin_sig[res_fibronectin_sig$log2FoldChange>0,] 
downregulated_fibronectin_sig <- res_fibronectin_sig[res_fibronectin_sig$log2FoldChange<0,] 
#write.csv(upregulated_fibronectin_sig, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/fibronectin_upregulated.csv", quote=FALSE)
#write.csv(downregulated_fibronectin_sig, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/fibronectin_downregulated.csv", quote=FALSE)


">>>>>>>>>>>>>>>>>>>>>>>Vitronectin>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

res_vitronectin <- results(dds, contrast = c("group", "vitronectin", "b27"))
res_vitronectin <- na.omit(res_vitronectin)
res_vitronectin_sig <- res_vitronectin[res_vitronectin$padj<0.05,]
res_vitronectin_sig <- res_vitronectin_sig[abs(res_vitronectin_sig$log2FoldChange)>0.58,]
upregulated_vitronectin_sig <- res_vitronectin_sig[res_vitronectin_sig$log2FoldChange>0,] 
downregulated_vitronectin_sig <- res_vitronectin_sig[res_vitronectin_sig$log2FoldChange<0,] 
#write.csv(upregulated_vitronectin_sig, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/vitronectin_upregulated.csv", quote=FALSE)
#write.csv(downregulated_vitronectin_sig, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/vitronectin_downregulated.csv", quote=FALSE)

">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
#getting genes that are upregulated in fibronectin, fibronectin and FBS



">>>>>>>>>>>>>>>>>>>>>Making Venn Diagrams>>>>>>>>>>>>>>>>>>>>>>>>>"

upregulated_venn <- venn.diagram(x <- list(row.names(upregulated_fbs_sig), row.names(upregulated_fibronectin_sig),
                       row.names(upregulated_vitronectin_sig)),
             category.names = c("FBS", "Fibronectin","Vitronectin"),
             filename = NULL,
             euler.d = TRUE,
             cat.fontfamily = "sans",
             scaled = TRUE,
             col="black",
             fontfamily = "sans",
             fill = c("#C0C0C0","#0066CC", "#FF9933"),
             #cat.col = c("#CCFFE5","#FFFFCC", "#FF99FF"),
             cat.cex = 1,
             margin = 0.1
)
ggsave(upregulated_venn, file = "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_2_venn.pdf",
       device = "pdf")


downregulated_venn <- venn.diagram(x <- list(row.names(downregulated_fbs_sig), row.names(downregulated_fibronectin_sig),
                       row.names(downregulated_vitronectin_sig)),
             category.names = c("FBS", "Fibronectin","Vitronectin"),
             filename = NULL,
             euler.d = TRUE,
             cat.fontfamily = "sans",
             scaled = TRUE,
             col="black",
             fontfamily = "sans",
             fill = c("#C0C0C0","#0066CC", "#FF9933"),
             #cat.col = c("#CCFFE5","#FFFFCC", "#FF99FF"),
             cat.cex = 1,
             margin = 0.1
)

ggsave(downregulated_venn, file = "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/downregulated_2_venn.pdf",
       device = "pdf")

">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Combining genes that are upregulated/downregulated in both ECMs of interest>>>>>>>>>>>>>>>>>"
fibronectin_downregulated_genes <- read.csv("/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/fibronectin_downregulated.csv")%>%
  select("X")
vitronectin_downregulated_genes <- read.csv("/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/vitronectin_downregulated.csv")%>%
  select("X")
fbs_downregulated_genes <- read.csv("/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/fbs_downregulated.csv")%>%
  select("X")
fibronectin_upregulated_genes <- read.csv("/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/fibronectin_upregulated.csv")%>%
  select("X")
vitronectin_upregulated_genes <- read.csv("/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/vitronectin_upregulated.csv")%>%
  select("X")
fbs_upregulated_genes <- read.csv("/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/fbs_upregulated.csv")%>%
  select("X")


#genes upregulated and downregulated in both fibronectin and vitronectin
upregulated_in_both_vitronectin_and_fibronectin <- merge(fibronectin_upregulated_genes, vitronectin_upregulated_genes)
downregulated_in_both_vitronectin_and_fibronectin <- merge(fibronectin_downregulated_genes,vitronectin_downregulated_genes)
upregulated_in_vitronectin_fibronectin_and_fbs <- merge(upregulated_in_both_vitronectin_and_fibronectin,fbs_upregulated_genes)
downregulated_in_vitronectin_fibronectin_and_fbs <- merge(downregulated_in_both_vitronectin_and_fibronectin,fbs_downregulated_genes)
upregulated_in_both_fbs_and_fibronectin <- merge(fibronectin_upregulated_genes, fbs_upregulated_genes)
downregulated_in_both_fbs_and_fibronectin <- merge(fibronectin_downregulated_genes,fbs_downregulated_genes)
upregulated_in_both_fbs_and_vitronectin <- merge(vitronectin_upregulated_genes, fbs_upregulated_genes)
downregulated_in_both_fbs_and_vitronectin <- merge(vitronectin_downregulated_genes,fbs_downregulated_genes)


#>>>>>>>>>>fibronectin only
upregulated_fbs_vitronectin_fibronectin_overlap <- data.frame(X =unique(c(upregulated_in_both_fbs_and_fibronectin$X, upregulated_in_both_vitronectin_and_fibronectin$X)))
upregulated_in_fibronectin_only <- anti_join(fibronectin_upregulated_genes,upregulated_fbs_vitronectin_fibronectin_overlap)
downregulated_fbs_vitronectin_fibronectin_overlap <- data.frame(X =unique(c(downregulated_in_both_fbs_and_fibronectin$X, downregulated_in_both_vitronectin_and_fibronectin$X)))
downregulated_in_fibronectin_only <- anti_join(fibronectin_downregulated_genes,downregulated_fbs_vitronectin_fibronectin_overlap)

#vitronectin only
upregulated_fbs_vitronectin_fibronectin_overlap_2 <- data.frame(X =unique(c(upregulated_in_both_fbs_and_vitronectin$X, upregulated_in_both_vitronectin_and_fibronectin$X)))
upregulated_in_vitronectin_only <- anti_join(vitronectin_upregulated_genes,upregulated_fbs_vitronectin_fibronectin_overlap_2)
downregulated_fbs_vitronectin_fibronectin_overlap_2 <- data.frame(X =unique(c(downregulated_in_both_fbs_and_vitronectin$X, downregulated_in_both_vitronectin_and_fibronectin$X)))
downregulated_in_vitronectin_only <- anti_join(vitronectin_downregulated_genes,downregulated_fbs_vitronectin_fibronectin_overlap_2)




#write.csv(upregulated_in_both_vitronectin_and_fibronectin,"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/upregulated_in_both_vitronectin_and_fibronectin.csv",quote=FALSE)
#write.csv(downregulated_in_both_vitronectin_and_fibronectin,"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/downregulated_in_both_vitronectin_and_fibronectin.csv",quote=FALSE)
#write.csv(downregulated_in_vitronectin_fibronectin_and_fbs,"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/downregulated_in_vitronectin_fibronectin_and_fbs",quote=FALSE)
#write.csv(upregulated_in_vitronectin_fibronectin_and_fbs,"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/upregulated_in_vitronectin_fibronectin_and_fbs",quote=FALSE)
#write.csv(upregulated_in_fibronectin_only, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/upregulated_in_fibronectin_only.csv",quote = F)
#write.csv(downregulated_in_fibronectin_only, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/downregulated_in_fibronectin_only.csv",quote = F)
#write.csv(upregulated_in_vitronectin_only, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/upregulated_in_vitronectin_only.csv",quote = F)
#write.csv(downregulated_in_vitronectin_only, "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/downregulated_in_vitronectin_only.csv",quote = F)



">>>>>>>>>>>>>>heatmap for the expression of integrin genes in fbs, fibronectin and vitronectin>>>>>>>>>>>>>>"
normalized_data <- vst(assay(dds))
n <- normalized_data[grep("ITG", rownames(normalized_data)),]
pheatmap(n, scale="row")


res_fbs <- results(dds, contrast = c("group", "fbs", "b27")) %>% 
  na.omit()
res_fbs_sig <- res_fbs[res_fbs$padj<0.05,]
res_fbs_int <- data.frame(res_fbs_sig[grep("ITG", rownames(res_fbs_sig)),"log2FoldChange", drop= FALSE])
res_fbs_int <- res_fbs_int %>% 
  mutate(int_genes = rownames(res_fbs_int), reg = log2FoldChange >0)
  
res_fibronectin <- results(dds, contrast = c("group", "fibronectin", "b27")) %>% 
  na.omit()
res_fibronectin_sig <- res_fibronectin[res_fibronectin$padj<0.05,]
res_fibronectin_int <- data.frame(res_fibronectin_sig[grep("ITG", rownames(res_fibronectin_sig)), "log2FoldChange",drop=FALSE])
res_fibronectin_int <- res_fibronectin_int %>% 
  mutate(int_genes = rownames(res_fibronectin_int), reg = log2FoldChange >0)


res_vitronectin <- results(dds, contrast = c("group", "vitronectin", "b27")) %>% 
  na.omit()
res_vitronectin_sig <- res_vitronectin[res_vitronectin$padj<0.05,]
res_vitronectin_int <- data.frame(res_vitronectin_sig[grep("ITG", rownames(res_vitronectin_sig)), "log2FoldChange", drop=FALSE])
res_vitronectin_int <- res_vitronectin_int %>% 
  mutate(int_genes = rownames(res_vitronectin_int), reg = log2FoldChange >0)

integrin_expression_plots <- function(data_, data_title){
  integrin_plot <- ggplot(data_, aes(y = fct_reorder(int_genes,log2FoldChange, .desc=TRUE),
                                     x = log2FoldChange, fill= log2FoldChange > 0))+ 
    geom_bar(stat="identity", width = 0.6) + 
    scale_fill_brewer(palette="Dark2") + 
    labs(y="Integrin genes", title = data_title, caption = "P-Adj value < 0.05")+
    theme_light()+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  return(integrin_plot)
}

integrin_expression_plots(res_vitronectin_int, "Vitronectin")


df_intregrin <- data.frame(all_intergins=unique(c("ITGB4", "ITGA3","ITGA5","ITGB1" ,"ITGB8","ITGAV","ITGA6","ITGA10","ITGB3BP",
                             "ITGB4","ITGA3","ITGA5","ITGB1","ITGB8","ITGAV","ITGA6","ITGB6","ITGB1BP1",
                             "ITGB4","ITGA3","ITGA5","ITGB1","ITGB8","ITGAV","ITGA6","ITGB1BP1","ITGA10","ITGB3BP" )))

colnames(res_fbs_int) <- c("log2FoldChange_fbs","int_genes","reg")
res_fbs_int_hs <- res_fbs_int[,c("log2FoldChange_fbs", "int_genes")]

colnames(res_vitronectin_int) <- c("log2FoldChange_vitronectin","int_genes","reg")
res_vitronectin_int_hs <- res_vitronectin_int[,c("log2FoldChange_vitronectin", "int_genes")]

colnames(res_fibronectin_int) <- c("log2FoldChange_fibronectin","int_genes","reg")
res_fibronectin_int_hs <- res_fibronectin_int[,c("log2FoldChange_fibronectin", "int_genes")]


int_merge_for_heatmap <- left_join(x = df_intregrin, y = res_fibronectin_int_hs, by = c("all_intergins" = "int_genes")) %>%
  left_join(y = res_fbs_int_hs, by = c("all_intergins"="int_genes"))%>% 
  left_join(y = res_vitronectin_int_hs, by = c("all_intergins"="int_genes"))


row.names(int_merge_for_heatmap) <- int_merge_for_heatmap$all_intergins
int_merge_for_heatmap <- int_merge_for_heatmap[, c("log2FoldChange_fibronectin","log2FoldChange_fbs",
                                                   "log2FoldChange_vitronectin")]
colnames(int_merge_for_heatmap) <- c("FBS", "Fibronectin", "Vitronectin")



int_merge_for_heatmap <- int_merge_for_heatmap[complete.cases(int_merge_for_heatmap),]
pheatmap(int_merge_for_heatmap, scale="row")



#write.csv(upregulated_in_both_vitronectin_and_fibronectin,"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/upregulated_in_both_vitronectin_and_fibronectin.csv",quote=FALSE)
#write.csv(downregulated_in_both_vitronectin_and_fibronectin,"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/downregulated_in_both_vitronectin_and_fibronectin.csv",quote=FALSE)

">>>>>>>>>>>>>>>>>>>>>>>>>>>>making the venn diagram for only DE genes in Fibronectin and vitronectin \
as justed by JI>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

upregulated_venn <- venn.diagram(x <- list(row.names(upregulated_fibronectin_sig),
                                           row.names(upregulated_vitronectin_sig)),
                                 category.names = c("Fibronectin","Vitronectin"),
                                 filename = NULL,
                                 euler.d = TRUE,
                                 cat.fontfamily = "sans",
                                 scaled = TRUE,
                                 col="black",
                                 fontfamily = "sans",
                                 fill = c("#0066CC", "#FF9933"),
                                 #cat.col = c("#CCFFE5","#FFFFCC", "#FF99FF"),
                                 cat.cex = 1,
                                 margin = 0.1
)
ggsave(upregulated_venn, file = "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_vitro_fibro_venn.pdf",
       device = "pdf")


downregulated_venn <- venn.diagram(x <- list(row.names(downregulated_fibronectin_sig),
                                             row.names(downregulated_vitronectin_sig)),
                                   category.names = c("Fibronectin","Vitronectin"),
                                   filename = NULL,
                                   euler.d = TRUE,
                                   cat.fontfamily = "sans",
                                   scaled = TRUE,
                                   col="black",
                                   fontfamily = "sans",
                                   fill = c("#0066CC", "#FF9933"),
                                   #cat.col = c("#CCFFE5","#FFFFCC", "#FF99FF"),
                                   cat.cex = 1,
                                   margin = 0.1
)

ggsave(downregulated_venn, file = "/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/downregulated_vitro_fibro_venn.pdf",
       device = "pdf")



">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>getting genes that were upregulated and downregulated only in fibronectin/vitronectin>>>>>>>>>>"
"
After makking venn diagram for genes diffrentially expressed in vitronectin and fibronectin \
I plan to get the list of genes that is uniquely upregulated/downregulated in  any of the two \
ECM. after this, the pathway analysis of this genes unique gene will be done. 

"

unique_ <- function(big_gene_list, intersection) {
  gene_list <- anti_join(big_gene_list, intersection, by = c("X"="X"))
  return(gene_list)
}

#write.csv(unique_(fibronectin_upregulated_genes,upregulated_in_both_vitronectin_and_fibronectin),"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/fibronectin_vitrronectin/upregulated_only_in_fibronectin.csv",quote= FALSE)
#write.csv(unique_(fibronectin_downregulated_genes,downregulated_in_both_vitronectin_and_fibronectin),"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/fibronectin_vitrronectin/downregulated_only_in_fibronectin.csv",quote= FALSE)
#write.csv(unique_(vitronectin_upregulated_genes,upregulated_in_both_vitronectin_and_fibronectin),"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/fibronectin_vitrronectin/upregulated_only_vitronectin.csv",quote= FALSE)
#write.csv(unique_(vitronectin_downregulated_genes,downregulated_in_both_vitronectin_and_fibronectin),"/Users/myname/Desktop/miapaca on ECM Project (Genomics)/03282020_deseq_analysis_after_combining_technical_replicates_using_the_result_funtion_without_the_lfcthreshold/upregulated_downregulated_combinations/fibronectin_vitrronectin/downregulated_only_in_vitronectin.csv",quote= FALSE)








