importData <- function(i) {
	# Read expression and phenotype values for each cancer type, i, and save to list

	exp_dir <- "/home/t.cri.cczysz/tcga/expression"
	phen_dir <- "/home/t.cri.cczysz/tcga/phen"
	files <- list.files(paste(exp_dir, i, sep='/'), pattern='*.data.txt')

	# Ignore first two header rows and import expression
	x <- read.table(paste(exp_dir, i, files[1],sep='/'), header=T, row.names=1, skip=2, sep='\t')
	# Import header only
	y <- read.table(paste(exp_dir, i, files[1],sep='/'), header=T, row.names=1, skip=0, sep='\t', nrows=2)
	colnames(x) <- colnames(y)
	rm(y)
	
	files <- list.files(paste(phen_dir, i, sep='/'), pattern='*.picked.txt')
	phen <- read.table(paste(phen_dir,i,files[1],sep='/'), header=F, row.names=1, skip=0, sep="\t")

	exp_ids <- tolower(apply(matrix(colnames(x), ncol=1), 1, function(x) {paste(unlist(strsplit(x, split='[.]'))[1:3], collapse='-')}))
	phen_ids <- as.character(t(phen)[,1])

	# Remove entries with duplicates
	colnames(x) <- exp_ids
	unique_exp <- x[, !(exp_ids%in%exp_ids[duplicated(exp_ids)])]

	phens <- t(phen)[phen[1,]%in%colnames(unique_exp),]
	exprs <- unique_exp[, colnames(unique_exp)%in%phens[,1]]
	fit <- performDE(exprs, phens)
	#return(list(phens, exprs))
	return(fit)

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

performDE <- function(expr, phen, filter=T, sv=F) {

	# Filtering:
	## RPKM >0.01 in 5% of samples

	expr <- as.matrix(expr)
	nsamples <- ncol(expr)
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

library(edgeR)
library(limma)
library(biomaRt)
library(ggplot2)
#library(goseq)
library(sva)


setwd('/home/t.cri.cczysz/gdc/')

load('exprs.Robj')
load('phens.Robj')

if (F) {
exp_dir <- "/home/t.cri.cczysz/tcga/expression"
phen_dir <- "/home/t.cri.cczysz/tcga/phen"

#cancer <- c("ACC")
cancer <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

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

de_outfile <- '/home/t.cri.cczysz/tcga/de_results.Robj'
#de_outfile <- '/home/t.cri.cczysz/tcga/de_results_sva.Robj'

if (!file.exists(de_outfile)) {
	#exprs <- list()
	#phens <- list()

	#out <- lapply(cancer, importData)
	et_list <- lapply(cancer, importData)
	names(et_list) <- cancer
	#names(out) <- cancer

	#exprs <- lapply(out, function(x) { return(x[[2]]) })
	#names(exprs) <- cancer
	#phens <- lapply(out, function(x) { return(x[[1]]) })
	#names(phens) <- cancer
	#rm(out)

	#et_list <- lapply(cancer, function(x, exprs, phens) {performDE(exprs[[x]], phens[[x]], T, F)})
	#names(et_list) <- cancer
	if (F) {
	for (i in cancer) {
		out <- importData(i)
		phens[[i]] <- out[[1]]
		exprs[[i]] <- out[[2]]
	}

	et_list <- list()
	for (i in cancer) {
		et <- performDE(exprs[[i]], phens[[i]], T, F)
		et_list[[i]] <- et
	}
	}

	#results <- lapply(et_list, topTable, number=Inf)
	save(et_list, file=de_outfile)

} else {load(file=de_outfile)}

# et_list[[cancer]] gives lmFit objects from lm of gene ~ sex test

# Perform GO enrichment analysis
#go.data <- list()

go_outfile <- '/home/t.cri.cczysz/tcga/go.Robj'
#if (!(file.exists(go_outfile)) | !(file.exists('/home/t.cri.cczysz/go_sig.Robj'))) {
if (!file.exists(go_outfile)) { 
	go.data <- performGO(et_list)
	#go.sig.data <- lapply(go.data, function(x) {subset(x,q.value<0.25)})

	save(go.data, file=go_outfile)
	#save(go.sig.data, file='/home/t.cri.cczysz/go_sig.Robj') 
} else {load(file=go_outfile)}

ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
cancer_summaries <- c('cancer', 'name', 'sample_size', 'nmale', 'nfemale', 'percentfemale', 'nsva','ngenes', 'nsiggenes', 'malebiased', 'percentmale', 'femalebiased', 'percentfemale')

summary_outfile = '/home/t.cri.cczysz/tcga/de_summary.Robj'

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

pdf(file='sig_genes.pdf', width=11, height=8)
	g <- ggplot(data=plot.data, aes(x=rownames(plot.data), y=as.numeric(nsiggenes)))
	g + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	g <- ggplot(data=plot.data, aes(x=rownames(plot.data), y=log10(as.numeric(nsiggenes))))
	g + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

if (F) {
for (cancer_type in names(out_list)) {
	f.out <- paste(cancer_type, 'results.csv', sep='.')
	write.csv(out_list[[cancer_type]], file=f.out, sep=',', row.names=T, col.names=T, quote=F)
}
}
}
