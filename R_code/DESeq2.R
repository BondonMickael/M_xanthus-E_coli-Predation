##----Import Library----
library(DESeq2)
library(apeglm)
library(gprofiler2)
library(enrichR)
library(ggplot2)

##----Import Data----
Predation = read.csv2(file = "~/Documents/Data/WT_vs_4preys_final.csv", 
                      header = TRUE, sep = ',', dec = '.', row.names="X")
Predation = na.omit(Predation)
Predation = subset(Predation, select = -c(K1, K2, K3, KC1, KC2, KC3, Y1, Y2, Y3)) # remove KC and K
apply(Predation, FUN = sum, MARGIN = 2)

##----Create coldata----
coldata = data.frame(colnames(Predation))
coldata$condition = coldata$colnames.Predation.
coldata$prey = coldata$colnames.Predation.

for(i in coldata$colnames.Predation.){
  if (grepl("Ec", i, fixed = TRUE)) {
    coldata$condition[coldata$colnames.Predation. == i] = "Predation"
    coldata$prey[coldata$colnames.Predation. == i] = "E.coli 1/3"
  }
  if (grepl("B", i, fixed = TRUE)) {
    coldata$condition[coldata$colnames.Predation. == i] = "Predation"
    coldata$prey[coldata$colnames.Predation. == i] = "B.subtilis"
  }
  if (grepl("C", i, fixed = TRUE)) {
    coldata$condition[coldata$colnames.Predation. == i] = "Predation"
    coldata$prey[coldata$colnames.Predation. == i] = "Caulobacter"
  }
  if (grepl("W", i, fixed = TRUE)) {
    coldata$condition[coldata$colnames.Predation. == i] = "Alone"
    coldata$prey[coldata$colnames.Predation. == i] = "ANone"
  }
  if (grepl("WC", i, fixed = TRUE)) {
    coldata$condition[coldata$colnames.Predation. == i] = "Predation"
    coldata$prey[coldata$colnames.Predation. == i] = "E.coli 1/1"
  }
}

coldata_final = coldata[,-1]
rownames(coldata_final) = coldata[,1]

##----Create/Perform DESeq2 (Diff Express Analysis)---
dds = DESeqDataSetFromMatrix(countData = Predation,
                             colData = coldata_final,
                             design= ~prey)
DEA = DESeq(dds)
resultsNames(DEA)
res_B = results(DEA, name = "prey_B.subtilis_vs_ANone")
res_C = results(DEA, name = "prey_Caulobacter_vs_ANone")
res_E = results(DEA, name = "prey_E.coli.1.1_vs_ANone")
res_E3 = results(DEA, name = "prey_E.coli.1.3_vs_ANone")
resLFC_B = lfcShrink(DEA, coef="prey_B.subtilis_vs_ANone", type="apeglm")
resLFC_C = lfcShrink(DEA, coef="prey_Caulobacter_vs_ANone", type="apeglm")
resLFC_E = lfcShrink(DEA, coef="prey_E.coli.1.1_vs_ANone", type="apeglm")
resLFC_E3 = lfcShrink(DEA, coef="prey_E.coli.1.3_vs_ANone", type="apeglm")

plotMA(res_B, ylim=c(-2,2), main = "B.subtilis")
plotMA(resLFC_B, ylim=c(-2,2), main = "B.subtilis shrink")

plotMA(res_C, ylim=c(-2,2), main = "Caulobacter")
plotMA(resLFC_C, ylim=c(-2,2), main = "Caulobacter shrink")

plotMA(res_E, ylim=c(-2,2), main = "E. coli")
plotMA(resLFC_E, ylim=c(-2,2), main = "E. coli shrink")

plotMA(res_E3, ylim=c(-2,2), main = "E. coli 1/3")
plotMA(resLFC_E3, ylim=c(-2,2), main = "E. coli 1/3 shrink")

#res[order(res$log2FoldChange,decreasing = TRUE),]

##----Volcano plot----
topT = as.data.frame(res_C)

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main= substitute(paste("Volcano plot for ", italic('M. xanthus'), " alone vs predation ", italic('Caulobacter'))), cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))
with(subset(topT, padj<0.05 | abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="cyan", cex=0.5))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
with(topT[rownames(topT) %in% genes_list, ], text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>4), cex=0.8, pos=3))
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>4), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>4), cex=0.8, pos=3))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)
#abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

##----Take Significantly Expressed Genes----
sum(na.omit(res_C$padj < 0.05 & abs(res_C$log2FoldChange) > 2))
Genes_id = res_C$padj < 0.05 & abs(res_C$log2FoldChange) > 2
Genes_C = na.omit(res_C@rownames[Genes_id]) 

sum(na.omit(res_B$padj < 0.05 & abs(res_B$log2FoldChange) > 2))
Genes_id = res_B$padj < 0.05 & abs(res_B$log2FoldChange) > 2
Genes_B = na.omit(res_B@rownames[Genes_id]) 

#sum(na.omit(res_E$padj < 0.05 & abs(res_E$log2FoldChange) > 2))
Genes_id = res_E$padj < 0.05 & abs(res_E$log2FoldChange) > 2
Genes_E = na.omit(res_E@rownames[Genes_id]) 

sum(na.omit(res_E3$padj < 0.05 & abs(res_E3$log2FoldChange) > 2))
Genes_id = res_E3$padj < 0.05 & abs(res_E3$log2FoldChange) > 2
Genes_E3 = na.omit(res_E3@rownames[Genes_id]) 

Liste_Genes = Genes_B[Genes_B %in% Genes_E3]
Liste_Genes = Liste_Genes[Liste_Genes %in% Genes_C]
lapply(Liste_Genes, write, "~/Documents/R/Liste_Genes", append=TRUE, ncolumns=1000)

##----Take only up or down----
#----Up----
Genes_id_up = res_C$padj < 0.05 & res_C$log2FoldChange > 2
Genes_C_up = na.omit(res_C@rownames[Genes_id_up]) 

Genes_id_up = res_B$padj < 0.05 & res_B$log2FoldChange > 2
Genes_B_up = na.omit(res_B@rownames[Genes_id_up]) 

Genes_id_up = res_E$padj < 0.05 & res_E$log2FoldChange > 2
Genes_E_up = na.omit(res_E@rownames[Genes_id_up]) 

Genes_id_up = res_E3$padj < 0.05 & res_E3$log2FoldChange > 2
Genes_E3_up = na.omit(res_E3@rownames[Genes_id_up]) 

Liste_Genes_up = Genes_B_up[Genes_B_up %in% Genes_E3_up]
Liste_Genes_up = Liste_Genes_up[Liste_Genes_up %in% Genes_C_up]
lapply(Liste_Genes_up, write, "~/Documents/R/Liste_Genes_up", append=TRUE, ncolumns=1000)

#----Down----
Genes_id_down = res_C$padj < 0.05 & res_C$log2FoldChange < -2
Genes_C_down = na.omit(res_C@rownames[Genes_id_down]) 

Genes_id_down = res_B$padj < 0.05 & res_B$log2FoldChange < -2
Genes_B_down = na.omit(res_B@rownames[Genes_id_down]) 

Genes_id_down = res_E$padj < 0.05 & res_E$log2FoldChange < -2
Genes_E_down = na.omit(res_E@rownames[Genes_id_down]) 

Genes_id_down = res_E3$padj < 0.05 & res_E3$log2FoldChange < -2
Genes_E3_down = na.omit(res_E3@rownames[Genes_id_down]) 

Liste_Genes_down = Genes_B_down[Genes_B_down %in% Genes_E3_down]
Liste_Genes_down = Genes_E3_down[Genes_E3_down %in% Genes_C_down]

#----Take list----
lapply(Genes_E_up, write, "~/Documents/R/Genes_List/Liste_Genes_E1_up", append=TRUE, ncolumns=1000)
lapply(Genes_E_down, write, "~/Documents/R/Genes_List/Liste_Genes_E1_down", append=TRUE, ncolumns=1000)

