# Gabriel Hoffman
# May 27, 2020
# 
# Template for gene set analyis

# Gene set testing
##################

library(data.table)
library(synapser)
synLogin()

# Read GMT file to list
# 	reads either .gmt or gmt.gz
read.gmt = function(file){
    if( ! grepl("(\\.gmt$)|(\\.gmt\\.gz$)", file)[1] ){
        stop("Pathway information must be a .gmt file")
    }
    geneSetDB = readLines(file)
    geneSetDB = strsplit(geneSetDB, "\t")
    names(geneSetDB) = sapply(geneSetDB, "[", 1)
    geneSetDB = lapply(geneSetDB, "[", -1:-2)
    geneSetDB = lapply(geneSetDB, function(x) {
        x[which(x != "")]
    })
    return(geneSetDB)
}

# Map from HGNC to ENSEMBL gene ids
convert_gene_list_to_ensembl = function( x, geneInfo){
	geneInfo[ x,]$ensembl_gene_id
}

# For ENSEMBL id ENSG00000279457.4, return ENSG00000279457
trim_ensembl_ids = function(x){
	gsub("(.*)\\.(.*)", "\\1", x)	
}

# Download GMT files storing gene sets
# combine multiple GMT files into a single list
syn_gs = c('syn22097899', 'syn22097897', 'syn22097895', 'syn22097894', 'syn22097898')

geneSets = lapply( syn_gs, function(id){
	read.gmt(synGet(id)$path)
})
geneSets = do.call(c, geneSets)

# get ENSEMBL to HGNC gene mapping
geneInfo = fread(synGet('syn21907981')$path)
setkey( geneInfo, hgnc_symbol)

# Convert gene sets from HGNC to ENSEMBL ids
geneSets_ens = lapply( geneSets, function(gs){
	convert_gene_list_to_ensembl( gs, geneInfo)
	})
names(geneSets_ens) = names(geneSets)


# Process RNA-seq Counts
#
# I'm starting from https://github.com/CommonMindConsortium/covarr-de/blob/9ba620c46d4082cf4589d65b576e542e85880b24/individual-cohorts/MSSM-Penn-Pitt_DLPFC_DE.Rmd#L139
# cd /hpc/users/hoffmg01/work/covarr-de/individual-cohorts
# library(pinnacle)
# rmarkdown::render('MSSM-Penn-Pitt_ACC_DE.Rmd')
# Note that the limma::camera() method for gene set analysis does
# 	not accept random effects

# Estimate Precision weights with voom
rownames(NEW.COUNTS) = gsub("\\.(\\d+)", "", rownames(NEW.COUNTS))
geneExpr = DGEList(NEW.COUNTS)
geneExpr = edgeR::calcNormFactors(geneExpr)

form = ~ Dx + ageOfDeath + RIN2 + RIN + IntronicRate + IntragenicRate + Institution + Reported_Gender
dsgn = model.matrix(form, COVARIATES)
vobj = voom(geneExpr, dsgn)
fit = lmFit( vobj, dsgn)


# map ENSEMBL ids in rows of vobj to gene sets
index = ids2indices( geneSets_ens, trim_ensembl_ids(rownames(vobj)))

# filter out gene sets with fewer than 20 entries
pass_filter = sapply(index, function(x) length(x) >= 20)
table(pass_filter)

# Perform gene set analysis with camera
# test coefficient: DxSCZ
# This method models the correlation structure within each gene set
# 	since the test statistics are not independent
# Setting inter.gene.cor=NA, estimates the mean correlation for each gene set
resGSA = camera( vobj, index[pass_filter], dsgn, contrast='DxSCZ', inter.gene.cor=NA, use.ranks=TRUE)


# recodeToSparseMatrix
library(Matrix)
idx_sub = index[pass_filter]#[1:3354]

n_genes = max(sapply(idx_sub, max))
M = matrix(0, length(idx_sub), n_genes )

for( i in 1:length(idx_sub) ){
	M[i,idx_sub[[i]]] = 1
}

M = M[,apply(M, 2, var)> 0]

library(poolr)
# meff(cor(t(M)), method="gao")
n = ncol(M)
A = scale(t(M)) / sqrt(n-1)
# meff(eigen=eigen(crossprod(A))$values, method="gao")

evs = svd(A)$d^2
meff(eigen=evs, method="gao")



fit = lmFit( vobj, dsgn)
fit = eBayes(fit)

table(topTable(fit, coef="DxSCZ", number=Inf)$adj.P.Val < 0.05)

tstat = topTable(fit, coef="DxSCZ", number=Inf, sort.by='none')$t

resGSA.pr = cameraPR( tstat, index[pass_filter], use.ranks=TRUE, inter.gene.cor=.01)


resGSA.pr2 = cameraPR( tstat^2, index[pass_filter], use.ranks=TRUE, inter.gene.cor=.01)




# Use GSEABase
##############

library(GSEABase)
# library(data.table)
library(pinnacle)
library(synapser)
synLogin()

syn_gs = c('syn22097899', 'syn22097897', 'syn22097895', 'syn22097894', 'syn22097898')

files = sapply( syn_gs, function(id){
    synGet(id)$path
    })
names(files) = gsub(".gmt(.gz+)", "", basename(files))


# Read GMT files into GeneSetCollection
# The gene identifiers are symbols so specify SymbolIdentifier()
gscObj = readGMT( files, SymbolIdentifier() )

saveRDS(gscObj, file="gscObj.RDS")

# convert Symbols to ENSEMBL
gscObj.ens = mapIdentifiers(gscObj, ENSEMBLIdentifier("org.Hs.eg"))


gsc.dt = recodeToDataTable( gscObj ) 

gsc.lst = recodeToList( gscObj )



# Then convert to sparseMatrix
gsc.sparse = recodeToSparseMatrix( gscObj )


C = getJaccardMatrix( gscObj[1:100] )

C = getJaccardMatrix( gsc.sparse[,1:500])


filterBySetSize( filterByGenes( gscObj, "NOD2"), 1)

# geneProperties
################

geneProperties = getGeneLengthAndGCContent(rownames(vobj), "hg38", mode="org.db")
geneProperties = data.frame(gene_id = rownames(geneProperties), geneProperties, stringsAsFactors=FALSE)

geneProperties = merge(geneProperties, data.frame(AveExpr = rowMeans(vobj$E)), by="row.names")[,-1]


# pinnacle
###########


# only keep genes in gene sets?

library(pinnacle)

coef = 'Reported_GenderMale'
L = as.matrix(getContrast(vobj, form, COVARIATES, coef))
colnames(L) = coef

# get mapping from getGeneMappings to genes
df_peakGeneMap = getGeneMapping( fit, L, target)

# Evaluate decomposition of correlation structure
corrStruct = evalCorrelationDecomp( vobj, fit, L[,coef] )
 
gscObj.filter = filterBySetSize(gscObj, 20)
gscObj.filter2 = filterByGenes( gscObj.filter, df_peakGeneMap$targetGene)

resEnrich = pinnacleEnrich( vobj, fit, L, coef, corrStruct, df_peakGeneMap, gscObj.filter2, quantileTransform=FALSE)

# account for gene properties
resEnrich.gp = pinnacleEnrich( vobj, fit, L, coef, corrStruct, df_peakGeneMap, gscObj.filter2, geneProperties = geneProperties, quantileTransform=FALSE)


table( topTable(eBayes(fit), coef=coef, number=Inf)$adj.P.Val < 0.05)

pdf("pinnacle_sex.pdf", width=8, height=12)
plot_spread( resEnrich.gp, nsets=10)
plot_spread( resEnrich.gp[resEnrich.gp$FDR.spread < 0.05,][-c(1:3),], nsets=50)
dev.off()



# score each individual
#######################

gsids = unique(c(resEnrich.gp$Geneset[1:50], 
        resEnrich.gp$Geneset[order(resEnrich.gp$p.shift)][1:50]))

scoreIndivs = lapply( 1:nrow(COVARIATES), function(i){

    cat("\r",i,"       ")

    COVARIATES$ID = rep(0, nrow(COVARIATES))
    COVARIATES$ID[i] = 1

    form = ~ Dx + ageOfDeath + RIN2 + RIN + IntronicRate + IntragenicRate + Institution + Reported_Gender + ID
    dsgn = model.matrix(form, COVARIATES)
    fit = lmFit( vobj, dsgn)

    coef = "ID"
    L = as.matrix(getContrast(vobj, form, COVARIATES, coef))
    colnames(L) = coef

    pinnacleEnrich( vobj, fit, L, coef, corrStruct, df_peakGeneMap, gscObj.filter2[gsids], geneProperties = geneProperties, verbose=FALSE, quantileTransform=FALSE)
    })

# consolidate scores: individuals by genesets
key = 'NABA_PROTEOGLYCANS'
df_score = lapply(scoreIndivs, function(x){
    with(subset(x,Geneset==key), beta.spread/se_beta.spread)
    })
df_score = do.call(rbind, df_score)

hist(df_score)


# gs = filterByGeneSetName(gscObj.filter2, 'Gao_Large_Intestine_24W_C2_MKI67pos_Progenitor')

# plotPinnacleResult( vobj, fit, L, coef, corrStruct, gs, target, df_peakGeneMap)


gscObj.filter2.ens =  mapIdentifiers(gscObj.filter2, ENSEMBLIdentifier("org.Hs.eg")) 

gscObj.filter2.ens.list = recodeToList( gscObj.filter2.ens )


index = ids2indices( gscObj.filter2.ens.list, trim_ensembl_ids(rownames(vobj)))


resGSA = camera( vobj, index, dsgn, contrast='Reported_GenderMale', inter.gene.cor=NA, use.ranks=TRUE)

head(resGSA)




