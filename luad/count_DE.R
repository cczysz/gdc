performDE <- function(expr, phen) {

	# First, filter genes with 0 count in all individuals
	use <- (rowSums(expr) > 0)
	expr <- expr[use,]

	# Filter low gene counts
	nmale <- sum(phen$Sex == 'male')
	maleUse <- (rowSums(expr[, phen$Sex=='male'] >= 6) > 10 ) 
	nfemale <- sum(phen$Sex == 'female')
	femaleUse <- (rowSums(expr[, phen$Sex=='female'] >= 6) > 10) 

	# Filter out genes without counts of 6 in at least 10 males or females
	expr <- expr[(maleUse | femaleUse),]
	print(dim(expr))

	dds <- DESeqDataSetFromMatrix(countData = expr, colData = phen, design = ~ Sex)
	dds <- estimateSizeFactors(dds)
	dat <- counts(dds, normalized=TRUE)
	#dds <- DESeq(dds, parallel=T)
	#return(dds)
	# Estimate surrogate variables
	#dds <- estimateSizeFactors(dds)
	#print(head(dat))
	mod <- model.matrix(~phen$Sex, colData(dds))
	mod0 <- model.matrix(~1, colData(dds))
	svseq <- svaseq(dat, mod, mod0, n.sv=2)
	# svseq$sv: matrix of SVs
	# svseq$n.sv: number of SVs


	ddssva <- dds
	ddssva$SV1 <- svseq$sv[,1]
	ddssva$SV2 <- svseq$sv[,2]
	design(ddssva) <- ~ SV1 + SV2 + Sex	
	ddssva <- DESeq(ddssva, parallel=T)
	#res <- results(ddssva, parallel=T)
	return(ddssva)
	
}


if (!require(DESeq2)) {
	library(BiocInstaller)
	biocLite('DESeq2')
	library(DESeq2)
}
#library(biomaRt)
#library(ggplot2)
#library(goseq)
library(sva)

library(BiocParallel)
register(MulticoreParam(4))

setwd('/home/t.cri.cczysz/gdc/luad')

load('luad_exprs.Robj')
load('luad_phens.Robj')

exp_stats <- expr_df[(nrow(expr_df)-4):nrow(expr_df),]
exprs <- expr_df[1:(nrow(expr_df)-5),]
thca_deseq <- performDE(exprs, cond)
save(thca_deseq, file='thca_de.Robj')
thca_results <- results(thca_deseq, parallel=T)
thca_results <- thca_results[order(thca_results$padj),]
save(thca_results, file='thca_results.Robj')
