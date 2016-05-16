#conversions <- read.table(file='kirp_htseq.tsv', header=T)

phens <- read.table(file='phens.txt', sep='\t', header=T)

files_in <- read.table(file='files.txt', header=F, stringsAsFactors=F)
dup_ind <- read.table(file='duplicated.txt', header=F)
omit <- (files_in$V4%in%dup_ind$V1)

files_in <- files_in[!omit, ]

phens_and_expr <- intersect(phens$bcr_patient_uuid, files_in$V4)
phens <- phens[phens$bcr_patient_uuid%in%phens_and_expr,]
# files_in <- files_in[!(duplicated(files_in$V4)),]
# Recreate phenotype and expression objects
idx <- match(files_in[,4], phens[,2])
idx2 <- match(phens[,2], files_in[,4])
phens <- phens[idx,]
phens$File <- files_in$V1

df_list <- lapply(phens$File, read.table, header=F, row.names=1)

expr_df <- data.frame(row.names=rownames(df_list[[1]]))
for (i in seq(length(df_list))) {
	expr_df <- data.frame(expr_df, df_list[[i]])
}
rm(df_list)
colnames(expr_df) <- phens[,1]
cond <- data.frame(row.names=phens[,2], Sex=factor(phens$gender))

save(expr_df, file='exprs.Robj')
save(cond, file='phens.Robj')
