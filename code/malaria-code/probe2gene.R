if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(c('ArrayExpress', 'hgu133plus2.db'))
library(ArrayExpress)
library(hgu133plus2.db)
library(dplyr)
library(tidyr)
library(limma)

df <- read.table('../../output/malaria_dataset/data_mat.csv', sep=',', header=T)
Y <- df[1]
X <- df[2:ncol(df)]

for (i in 1:ncol(X)) {
    if (substr(colnames(X)[i], 1, 1) == "X") {
        colnames(X)[i] <- substr(colnames(X)[i], 2, nchar(colnames(X)[i]))
    }
}

anno <- AnnotationDbi::select(hgu133plus2.db,
                              keys = (colnames(X)),
                              columns = c("SYMBOL"), keytype = "PROBEID")

anno <- subset(anno, !is.na(SYMBOL))

anno_grouped <- group_by(anno, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

head(anno_summarized)

anno_kept <- filter(anno_summarized, no_of_matches == 1)

head(anno_kept)

ids_to_include <- (colnames(X) %in% anno_kept$PROBEID)

X_filt <- bind_cols(PROBEID = colnames(X)[ids_to_include], tibble(EXPR = t(X[,ids_to_include,drop=T])))

X_filt <- left_join(X_filt, anno)

X_filt <- select(X_filt, -PROBEID)

X_gene <- limma::avereps(X_filt$EXPR, ID = X_filt$SYMBOL)

X_gene <- t(X_gene)

data_gene <- cbind(Y, X_gene)

write.table(data_gene, '../../output/malaria_dataset/data_gene_mat.csv', row.names=FALSE, sep=',')
