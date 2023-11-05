
#module load R/4.2.0

#cd /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/featurecounts

R

counts <- read.table("gene_name_counts_s2_rf.txt",header=TRUE) 
x<- colnames(counts)
x <- gsub("X.projects.b1042.Shukla_lab.Daniela.projects.11.RNA.seq_Cg1_invitro.star.", "", x)
x <- gsub(".Aligned.out.sorted.PCR.bam", "", x)
colnames(counts) <- x
rownames(counts) <- as.vector(counts[,1])
counts <- counts[,6:ncol(counts)]

test_part1 = counts[,2:ncol(counts)]/(counts[,1]/1000)
test_part2 = (colSums(test_part1))/1000000
test_part3<-((t(apply(test_part1, 1, "/", test_part2))))
tpm <- test_part3

genenamee <- (paste(rownames(tpm), "'", sep=""))
tpm2 <- data.frame(genenamee, tpm)

write.table(tpm2, "20230710_gene_names_tpm_all_genes.csv", col.names=TRUE, row.names=TRUE, sep=",", quote=FALSE)

#
counts2 <- counts[,2:ncol(counts)]


epn.countTable = counts2
colnames(epn.countTable) = colnames(counts2)

conditions <- c("Ctrl","Ctrl","Ctrl","KO", "KO","KO")
conditions2 <- unique(conditions)
epn.Design = data.frame(
    row.names = colnames( epn.countTable ),
	condition = conditions
)


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = epn.Design,
                              design = ~ condition)


dds <- dds[ rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)

rld <- rlog(dds)

pdf("20230710_PCA.pdf")
plotPCA(rld)
dev.off()

v1 <- c("KO")
v2 <- c("Ctrl")

normcounts <- (counts(dds, normalized=TRUE))

normcountse <- (counts(dds, normalized=TRUE))



res <- NULL
for (i in 1:1){
res[i] <- list(results(dds,alpha=0.05,  contrast=c("condition",paste(v1[i]),paste(v2[i]))))
}

rescut <- NULL

for (i in 1:1){
res[[i]][(is.na(res[[i]][,6])),6]=2
rescut[i] <- list((res[[i]][abs(res[[i]][,2])>1 & res[[i]][,6] < 0.05,]))
}

tempall <- NULL
tempders <- NULL
for (i in 1:1){
tempall <- NULL
tempders <- NULL
tempders <- rescut[[i]]
tempall <- data.frame(normcountse,res[[i]])
write.table(tempall,paste("20230710_TE_All_regions_DESEQ_",v1[i],"_vs_",v2[i],".csv", sep="" ), row.names=TRUE, col.names=TRUE, quote=FALSE, sep=",")
tempders <- tempall[rownames(tempall) %in% rownames(rescut[[i]]),]
write.table(tempders,paste("20230710_TE_DEGs_",v1[i],"_vs_",v2[i],".csv", sep="" ), row.names=TRUE, col.names=TRUE, quote=FALSE, sep=",")

}

