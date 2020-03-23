# Gabriel Hoffman
# March 22, 2020

# Import GTF and save gene annotations as TSV file

library(rtracklayer)
library(synapser)
library(data.table)

synLogin()

# GENCODE 30 GTF used for RNA-seq quantification 
file = "/sc/hydra/projects/PBG/REFERENCES/GRCh38/Gencode/release_30/gencode.v30.primary_assembly.annotation.gtf"

GTF = import( file )

df_GTF = data.table(data.frame(GTF, stringsAsFactors=FALSE))

# Gene
#######

df_gene = df_GTF[type=='gene',]
setkey(df_gene, "gene_id")

ensGenes = unique( df_gene$gene_id )

res_gene = lapply( ensGenes, function(ensGene){
	df_gene[gene_id==ensGene,data.frame(seqnames, start, end, strand, type, gene_id, gene_name, gene_type, stringsAsFactors=FALSE)]
})
res_gene = data.table(do.call("rbind", res_gene))

res_gene[,TSS:=ifelse(strand=='+', start, end)]
res_gene[,TES:=ifelse(strand=='+', end, start)]

write.table(res_gene, file=gzfile("gencode_v30_gene_annotation.tsv.gz"), row.names=FALSE, sep="\t", quote=FALSE)

# Transcript
############

df_tr = df_GTF[type=='transcript',]
setkey(df_tr, "gene_id")

ensGenes = unique( df_tr$gene_id )

res_tr = lapply( ensGenes[1:1000], function(ensGene){
	df_tr[gene_id==ensGene,data.frame(seqnames, start, end, strand, type, gene_id, gene_name, gene_type, transcript_id, transcript_type, stringsAsFactors=FALSE)]
})
res_tr = do.call("rbind", res_tr)

res_tr[,TSS:=ifelse(strand=='+', start, end)]
res_tr[,TES:=ifelse(strand=='+', end, start)]

write.table(res_tr, file=gzfile("gencode_v30_transcript_annotation.tsv.gz"), row.names=FALSE, sep="\t", quote=FALSE)




# Upload Synapse results
########################

# gene
activity <- Activity(
 'GENCODE 30 gene annotation',
 description='GENCODE 30 gene annotation',
 used=c('syn21814507'),
 executed='https://github.com/GabrielHoffman/cmc_analysis/blob/master/src/get_gene_annotations.R')
file <- File('gencode_v30_gene_annotation.tsv.gz', description='GENCODE 30 gene annotation', parent='syn21814506')
file <- synStore(file, activity=activity)


# transcript
activity <- Activity(
 'GENCODE 30 transcript annotation',
 description='GENCODE 30 transcript annotation',
 used=c('syn21814507'),
 executed='https://github.com/GabrielHoffman/cmc_analysis/blob/master/src/get_gene_annotations.R')
file <- File('gencode_v30_transcript_annotation.tsv.gz', description='GENCODE 30 transcript annotation', parent='syn21814506')
file <- synStore(file, activity=activity)


# library(xlsx)

# df = read.xlsx("~/Downloads/GREX_PrediXcan_TWAS_SMR_GeneList_11.04.19.xlsx", 1)
# df$DISORDER = as.character(df$DISORDER)
# df$Gene = as.character(df$Gene)

# ids =  unique(gsub(' ', '', unlist(strsplit(ids, ','))))

# res_list = lapply(ids, function(id){
# 	df[grep(id,df$DISORDER),1]
# 	})
# names(res_list) = ids






# # GENCODE 30 Corresponds to ENSEMBL v99

# library(EnsDb.Hsapiens.v96)

# ensdb = EnsDb.Hsapiens.v86

# columns(ensdb)

# listColumns(ensdb)


# select(ensdb, keys = list(GeneNameFilter(c("BCL2", "BCL2L11"))),
# 	columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE", "TXSEQSTART",
#                        "TXSEQEND", "SEQNAME", "SEQSTRAND", "GENESEQSTART", "GENESEQEND"))

