importData <- function(i) {
	# Read expression and phenotype values for each cancer type, i, and save to list
	print(i)
	exp_dir <- "/home/t.cri.cczysz/tcga/mirna_expression"
	phen_dir <- "/home/t.cri.cczysz/tcga/phen"

	files <- list.files(paste(exp_dir, i, sep='/'), pattern='*.data.txt')

	# Ignore first two header rows and import expression
	# First column gives miRNA_ID
	# Each individual sample has 3 columns: read_count, reads_per_million, cross-mapped
	x <- read.table(paste(exp_dir, i, files[1],sep='/'), header=T, row.names=1, skip=2, sep='\t')
	# Import header only
	y <- read.table(paste(exp_dir, i, files[1],sep='/'), header=F, row.names=1, skip=0, sep='\t', nrows=1)
	y <- unique(as.character(y))

	rpm <- data.frame(x[, seq(2, ncol(x), 3)])
	colnames(rpm) <- y
	counts <- data.frame(x[, seq(1, ncol(x), 3)])
	colnames(counts) <- y

	files <- list.files(paste(phen_dir, i, sep='/'), pattern='*.picked.txt')
	phen <- read.table(paste(phen_dir,i,files[1],sep='/'), header=F, row.names=1, skip=0, sep="\t")
	exp_ids <- tolower(apply(matrix(colnames(counts), ncol=1), 1, function(x) {paste(unlist(strsplit(x, split='[-]'))[1:3], collapse='-')}))
	colnames(counts) <- exp_ids
	
	phen_ids <- as.character(t(phen)[,1])
	unique_exp <- counts[, !(exp_ids%in%exp_ids[duplicated(exp_ids)])]

	phens <- t(phen)[phen[1,]%in%colnames(unique_exp),]
	exprs <- unique_exp[, colnames(unique_exp)%in%phens[,1]]

	phen <- data.frame(row.names = phens[,1], Sex = factor(phens[,'gender']))
	# Reorder phen data frame by order of expression matrix
	phen <- phen[match(colnames(exprs), rownames(phen)),, drop=F]
	#return(list(exprs, phen))
	if (length(levels(phen$Sex)) == 2) {
		fit <- performDE(exprs, phen) 
	} else { fit <- NULL }
	return(fit)
	#return(list(exprs, phen))

}

performDE <- function(expr, phen) {

	# First, filter genes with 0 count in all individuals
	use <- (rowSums(expr) > 0)
	expr <- expr[use,]

	# Filter low gene counts
	nmale <- sum(phen$Sex == 'male')
	maleUse <- (rowSums(expr[, phen$Sex=='male'] > 10) >= (0.66 * nmale))
	nfemale <- sum(phen$Sex == 'female')
	femaleUse <- (rowSums(expr[, phen$Sex=='female'] > 10) >= (0.66 * nfemale))

	expr <- expr[(maleUse & femaleUse),]
	print(dim(expr))
	if (F) {
	gt10 <- expr > 10 # True if gene x individual count greater than 10
	rs <- rowSums(gt10) # For each gene, number of individuals w/ count > 10
	use <- (rs >= (0.95*ncol(expr))) # Test genes with 95% individuals having counts > 10 
	expr <- expr[use,]
	print(dim(expr))
	}

	dds <- DESeqDataSetFromMatrix(countData = expr, colData = phen, design = ~ Sex)
	dds <- DESeq(dds, parallel=T)
	return(dds)
	if (F) {	
	# Estimate surrogate variables
	#dds <- estimateSizeFactors(dds)
	dat <- counts(dds, normalized=TRUE)
	#print(head(dat))
	mod <- model.matrix(~phen$Sex, colData(dds))
	mod0 <- model.matrix(~1, colData(dds))
	svseq <- svaseq(dat, mod, mod0, n.sv=2)
	# svseq$sv: matrix of SVs
	# svseq$n.sv: number of SVs


	ddssva <- dds
	ddssva$SV1 <- svseq$sv[,1]
	ddssva$SV2 <- svseq$sv[,2]
	design(ddssva) <- ~ SV1 + Sex	
	ddssva <- DESeq(ddssva)
	res <- results(ddssva)
	return(res)
	}
	#condition <- as.factor(phen$Sex)
	#cds <- newCountDataSet(expr, condition)	
	#cds <- estimateSizeFactors(cds)

	#cds <- estimateDispersions(cds, method='blind', sharingMode='fit-only')
	#res <- nbinomTest(cds, 'male', 'female')
	#return(res)
if (F) {
	# Filtering:
	## RPKM >0.01 in 5% of samples

	#expr <- as.matrix(expr)
	#nsamples <- ncol(expr)
	#log2expr <- log2(expr+1)	

	# Filter lowly expressed genes
	expr <- expr[rowSums(expr)>0,] # Remove genes with 0 RPKM in all samples
	lowexp <- expr <= 0.1	
	lowexpr <- rowSums(lowexp) > (0.95 * nsamples)
	log2expr <- log2(expr[!lowexpr,]+1)

	sex <- as.numeric(phen[,'gender']=='male')

	if (sv & (sum(sex)==0 | sum(sex) == length(sex))) {
		sv <- F } 
	if (sv) {
		mod = model.matrix(~0 + as.factor(sex))
		mod0 = model.matrix(~0 + rep(1, nsamples))

		svobj <- sva(log2expr, mod, mod0)
		modSv <- cbind(mod, svobj$sv)

		design <- modSv
		contrast.matrix <- c(-1, 1, rep(0, svobj$n.sv))
		y <- DGEList(log2expr, remove.zeros=T)
		v <- voom(y, design)
		fit <- lmFit(v, design)
		contrast.fit <- contrasts.fit(fit, contrast.matrix)
		fit <- eBayes(contrast.fit)
		fit$nsamples <- ncol(log2expr)
	} else {
		design <- model.matrix(~sex)
		y <- DGEList(log2expr, remove.zeros=T)
		v <- voom(y, design)
		fit <- lmFit(v, design)
		fit <- eBayes(fit)
		#contrast.matrix <- c(-1, 1)
	}

	fit$nsamples <- ncol(log2expr)
	return(fit)
}
}

performGO <- function(et_list) {
	go.data <- list()
	for (i in names(et_list)) {
		tmp <- et_list[[i]]
		ttable <- topTable(tmp, number=Inf)
		sig <- topTable(tmp, number=Inf, p.value=0.05)
		all_genes <- rownames(topTable(tmp, number=Inf))
		go_in <- as.numeric(all_genes%in%rownames(sig))
		all_genes <- apply(matrix(all_genes, ncol=1), 1, function (x) {unlist(strsplit(x, '[|]'))[2]})
		names(go_in) <- all_genes
		pwf <- nullp(go_in, "hg19", "refGene")
		GO.wall <- goseq(pwf, "hg19", "refGene", method='Sampling', repcnt=2500)
		GO.wall$q.value <- p.adjust(GO.wall$over_represented_pvalue,method="fdr")
		go.data[[i]] <- GO.wall
	}
	return(go.data)
}

#library(DESeq2)
if (!require(DESeq2)) {
	library(BiocInstaller)
	biocLite('DESeq2')
	library(DESeq2)
}
library(biomaRt)
library(ggplot2)
#library(goseq)
library(sva)
library("BiocParallel")
register(MulticoreParam(4))

setwd('/home/t.cri.cczysz/gdc/')

load('exprs.Robj')
load('phens.Robj')

exp_stats <- expr_df[(nrow(expr_df)-4):nrow(expr_df),]
exprs <- expr_df[1:(nrow(expr_df)-5),]
x <- performDE(exprs, cond)
save(x, file='de.Robj')
y <- results(x, parallel=T)
save(y, file='results.Robj')

q()
if (F) {
exp_dir <- "/home/t.cri.cczysz/tcga/mirna_expression"
phen_dir <- "/home/t.cri.cczysz/tcga/phen"

#cancer <- c("ACC")
cancer <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

full_names <- c('Adrenocortical Carcinoma', 
	'Breast Lobular Carcinoma',
	'Breast Invasive Carcinoma',
	'Cervical Squamous Cell Carcinoma',
	'Cholangiocarcinoma',
	'Colorectal Adenocarcinoma', 
	'Diffuse Large B-cell',
	'Esophageal Carcinoma',
	'Glioblastoma Multiforme',
	'Head/Neck Squamous Cell Carcinoma',
	'Kidney Chromophobe',
	'Kidney Reneal Clear Cell',
	'Kidney Renal Papillary Cell Carcinoma',
	'Acute Myeloid Lukemia',
	'Lower Grade Glioma',
	'Liver Hepatocellular Carcinoma',
	'Lung Adenocarcinoma',
	'Lung Squaous Cell Carcinoma',
	'Ovarian',
	'Pancreatic Adenocarcinoma',
	'Pheochromocytoma and Paraganglioma',
	'Prostate Adenocarcinoma',
	'Rectum Adenocarcinoma',
	'Sarcoma',
	'Stomach Adenocarcinoma',
	'Testicular Germ Cell',
	'Thyroid Cancer',
	'Thymoma',
	'Uterine Corpus Endometrial Carcinoma',
	'Uterine Carcinosarcoma',
	'Uveal Melanoma'
)
exp_dir <- "/home/t.cri.cczysz/tcga/mirna_expression"
phen_dir <- "/home/t.cri.cczysz/tcga/phen"

de_outfile <- '/home/t.cri.cczysz/tcga/mirna_results/de_results.Robj'
#de_outfile <- '/home/t.cri.cczysz/tcga/de_results_sva.Robj'

deseq_list <- NULL
deseq_res_list <- NULL
#if (file.exists(de_outfile)) {
if (T) {
	#exprs <- list()
	#phens <- list()

	#out <- lapply(cancer, importData)
	deseq_list <- lapply(cancer, importData)
	names(deseq_list) <- cancer
	save(deseq_list, file=de_outfile)

} else {load(file=de_outfile)}

deseq_res_list <- lapply(deseq_list, results)
# et_list[[cancer]] gives lmFit objects from lm of gene ~ sex test

# Perform GO enrichment analysis
#go.data <- list()

pdf('mirna_ma.pdf', width=12, height=10)
for (i in seq(length(et_list))) {
	if (!is.null(et_list[[i]])) {plotMA(et_list[[i]], main=full_names[i], alpha=0.05)}
}
dev.off()
q()

if (F) {
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
cancer_summaries <- c('cancer', 'name', 'sample_size', 'nmale', 'nfemale', 'percentfemale', 'nsva','ngenes', 'nsiggenes', 'malebiased', 'percentmale', 'femalebiased', 'percentfemale')

summary_outfile = '/home/t.cri.cczysz/tcga/mirna_results/de_summary.Robj'

if (!file.exists(summary_outfile)) {
	out_list <- list()
	for (i in seq(length(et_list))) {
		cancer_type <- names(et_list)[i]
		nmale <- sum(et_list[[i]]$design[,2])
		nfemale <- sum(!et_list[[i]]$design[,2])
		nsva <- ncol(et_list[[i]]$design) - 2

		all_genes <- topTable(et_list[[i]], number=Inf)
		sig_genes <- topTable(et_list[[i]], number=Inf, p.value=0.05)

		ngenes <- nrow(all_genes)	
		nsig <- nrow(sig_genes)
		male_bias <- sum(sig_genes$logFC > 0)
		female_bias <- sum(sig_genes$logFC < 0)

		cancer_summaries <- rbind(cancer_summaries, c(cancer_type, full_names[i], nmale+nfemale, nmale, nfemale, nfemale / (nmale + nfemale), nsva, ngenes, nsig, male_bias, signif(100 * male_bias/nsig, 2), female_bias, signif(100*female_bias/nsig, 2)))

		sig_symbols <- apply(matrix(rownames(sig_genes), ncol=1), 1, function(x) {unlist(strsplit(x, split='[|]'))[1]})
		sig_entrez <- apply(matrix(rownames(sig_genes), ncol=1), 1, function(x) {unlist(strsplit(x, split='[|]'))[2]})

		if (length(sig_symbols)==0) {
			#out_list[[cancer_type]] <- NA
			next
		}

		biomart.results=getBM(ensembl, attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene", "chromosome_name"),filters="entrezgene",values=sig_entrez)

		# Reorder biomart results by entrez id
		x <- biomart.results[match(sig_entrez, biomart.results$entrezgene),]
		out <- cbind(x, logFC = signif(sig_genes$logFC,2), p.val = signif(sig_genes$P.Value,2), q.value = signif(sig_genes$adj.P.Val,2))
		rownames(out) <- rownames(sig_genes)
		out_list[[cancer_type]] <- out
	}       
	save(out_list, file=summary_outfile)
} else {load(file=summary_outfile)}
	
write.table(cancer_summaries, file='cancer_summaries.csv', sep=',', row.names=F, col.names=F, quote=F)

#plot.data <- data.frame(meta_list[-1,-1])
#colnames(plot.data) <- meta_list[1,][-1]
#rownames(plot.data) <- meta_list[,1][-1]

if (F) {
pdf(file='sig_genes.pdf', width=11, height=8)
	g <- ggplot(data=plot.data, aes(x=rownames(plot.data), y=as.numeric(nsiggenes)))
	g + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	g <- ggplot(data=plot.data, aes(x=rownames(plot.data), y=log10(as.numeric(nsiggenes))))
	g + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

for (cancer_type in names(out_list)) {
	f.out <- paste(cancer_type, 'results.csv', sep='.')
	write.csv(out_list[[cancer_type]], file=f.out, sep=',', row.names=T, col.names=T, quote=F)
}
}
}
}
