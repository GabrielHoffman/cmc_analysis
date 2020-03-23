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

res_gene = lapply( ensGenes[1:10], function(ensGene){
	df_gene[gene_id==ensGene,data.frame(seqnames, start, end, strand, type, gene_id, gene_name, gene_type)]
})
res_gene = do.call("rbind", res_gene)

write.table(res_gene, file=gzfile("gencode_v30_gene_annoation.tsv.gz"), row.names=FALSE, sep="\t", quote=FALSE)



# Transcript
############

df_tr = df_GTF[type=='transcript',]
setkey(df_tr, "gene_id")

ensGenes = unique( df_tr$gene_id )

res_tr = lapply( ensGenes[1:10], function(ensGene){
	df_tr[gene_id==ensGene,data.frame(seqnames, start, end, strand, type, gene_id, gene_name, gene_type, transcript_id, transcript_type)]
})
res_tr = do.call("rbind", res_tr)

write.table(res_tr, file=gzfile("gencode_v30_transcript_annoation.tsv.gz"), row.names=FALSE, sep="\t", quote=FALSE)




# Upload Synapse results
########################

 # Adding files with provenance:
 #
 # A synapse entity *syn1906480* contains data
 # entity *syn1917825* contains code
 #
 activity <- Activity(
     'GENCODE 30 gene annotation',
     description='GENCODE 30 gene annotation',
     used=c('syn1906480', 'http://data_r_us.com/fancy/data.txt'),
     executed='syn1917825')
 file <- File('gencode_v30_gene_annoation.tsv.gz', description='GENCODE 30 gene annotation', parent=project)
 file <- synStore(file, activity=activity)










# # GENCODE 30 Corresponds to ENSEMBL v99

# library(EnsDb.Hsapiens.v96)

# ensdb = EnsDb.Hsapiens.v86

# columns(ensdb)

# listColumns(ensdb)


# select(ensdb, keys = list(GeneNameFilter(c("BCL2", "BCL2L11"))),
# 	columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE", "TXSEQSTART",
#                        "TXSEQEND", "SEQNAME", "SEQSTRAND", "GENESEQSTART", "GENESEQEND"))

