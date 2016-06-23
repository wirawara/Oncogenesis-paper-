### R code from vignette source 'microRNA.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: load_packages
###################################################
library(aroma.affymetrix)
library(biomaRt)
library(limma)
library(gdata)
library(WriteXLS)



###################################################
### code chunk number 2: point_to_relevant_files
###################################################

HOME <- getwd()
setwd("/home/Shared/data/array/marra_exon/aroma/")

chipType <- "HuEx-1_0-st-v2"
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="coreR3,A20071112,EP")
cs <- AffymetrixCelSet$byName("marra", cdf=cdf)

setCdf(cs, cdf)



###################################################
### code chunk number 3: pre_processing
###################################################
bc <- RmaBackgroundCorrection(cs)
csBC <- process(bc,verbose=verbose)
qn <- QuantileNormalization(csBC, typesToUpdate="pm")
csN <- process(qn, verbose=verbose)



###################################################
### code chunk number 4: fit_transcript_probe_level_models
###################################################
plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
fit(plmTr, verbose=verbose)



###################################################
### code chunk number 5: quality_checks
###################################################
qamTr <- QualityAssessmentModel(plmTr)
pdf(file.path(HOME,"marra_quality.pdf"),w=11,h=8.5)
par(mai=c(2.5,0.82,0.82,0.42))
plotNuse(qamTr)


###################################################
### code chunk number 6: dev_off1
###################################################
dummy <- dev.off()


###################################################
### code chunk number 7: microRNA.Rnw:128-130
###################################################
pdf(file.path(HOME,"marra_rle_quality.pdf"),w=11,h=8.5)
par(mai=c(2.5,0.82,0.82,0.42))


###################################################
### code chunk number 8: microRNA.Rnw:133-134
###################################################
plotRle(qamTr)


###################################################
### code chunk number 9: dev_off2
###################################################
x <- dev.off()


###################################################
### code chunk number 10: quality_checks
###################################################
units <- readCdfUnits(getPathname(cdf))
units.probes <- lapply(units, function(x) names(x$groups))
n <- sapply(units.probes,length)

ids <- data.frame(transcript_cluster_id=rep(names(units.probes),n),
                  probeset_id=unlist(units.probes, use.names=FALSE),
                  stringsAsFactors=FALSE)



###################################################
### code chunk number 11: estimates
###################################################
cesTr <- getChipEffectSet(plmTr)
trFit <- extractDataFrame(cesTr, addNames=TRUE, verbose=verbose)
dm <- log2(as.matrix(trFit[,-c(1:5)]))
rownames(dm) <- trFit$groupName
dm<- cbind(dm,as.numeric(trFit$unitName))
colnames(dm)[85]="unitName"


###################################################
### code chunk number 12: biomart_fetching
###################################################
mart= useMart("ensembl")
mart = useDataset("hsapiens_gene_ensembl",mart = mart)
annotation = getBM(attributes = c("affy_huex_1_0_st_v2","hgnc_symbol","strand","transcript_start","chromosome_name"), 
	filters="affy_huex_1_0_st_v2", values = c(ids$probeset_id), mart=mart)
 


###################################################
### code chunk number 13: annotation
###################################################
o = order(annotation$hgnc_symbol, decreasing=TRUE)
annotation = annotation[o,]


m = match(annotation$affy_huex_1_0_st_v2, ids$probeset_id)
annotation$transcript_cluster_id = ids$transcript_cluster_id[m]
                 


###################################################
### code chunk number 14: loading_targets
###################################################
tarDir = "/home/andreak/marra_exon/miRNAtargets/"

mirs = c(195,296,497)
f = file.path(tarDir,paste("miRNAwalk",mirs,".xls",sep=""))
names(f) = paste("miR",mirs,sep="")

targets = lapply(f, function(u) {
	d = read.xls(u)
	m1 = match(d$gene_name, annotation$hgnc_symbol)
	b = cbind(d, annotation[m1,])
	na.omit(b[!duplicated(b$transcript_cluster_id),])
})



###################################################
### code chunk number 15: read_polyps_and_nonpolyps
###################################################
# splitting metadata
metadata = read.xls(paste(tarDir,"non+poly.xls",sep=""))
metadata$PatientID = factor(metadata$PatientID)
rownames(metadata) = gsub(".CEL  ","",metadata$CEL)

count = table(metadata$PatientID, metadata$TissueType)
splitting = apply(count, 2, function(x) rownames(count)[which(x>0)])
poly = metadata[metadata$PatientID %in% splitting$pedun_polyp,]
nonpoly = metadata[metadata$PatientID %in% splitting$nonpolyp_lesion,]
data = list(full = metadata, nonpoly = nonpoly, poly = poly)




###################################################
### code chunk number 16: making_gene_anno
###################################################
m = match(rownames(dm), annotation$transcript_cluster_id)
gAnno = annotation[m,c(6,2:5)]



###################################################
### code chunk number 17: modeling
###################################################

dm.tissue = storage = vector("list",3); 
names(dm.tissue) = names(storage) = names(data)
for(j in 1:length(data)){
  lesion = (data[[j]]$TissueType != "normal_mucosa")+0                      
  pat = factor(as.character(data[[j]]$PatientID))                           
  design = model.matrix(~ pat + lesion,data=data[[j]]) 
  m = match(rownames(data[[j]]), colnames(dm))
  dm.new = dm[,m]
  dm.tissue[[j]] = dm.new 	
  fit = lmFit(dm.new, design,method = "robust")
  fit$genes = gAnno
  fB = eBayes(fit)
  storage[[j]] = lapply(targets, function(u) {
	ids = u$transcript_cluster_id
	topTable(fB[ids,],coef="lesion",number=length(ids))
  })
}



###################################################
### code chunk number 18: writing_a_file
###################################################
st = unlist(storage, recursive=FALSE)	
z = rep(names(st),sapply(st,nrow))
zz = strsplit(z,"\\.")
end = cbind(dataset=sapply(zz,.subset,1),mir=sapply(zz,.subset,2),do.call(rbind,st))
df = split(end, end$mir)

WriteXLS("df", ExcelFileName="miRNAs.xls", SheetNames=c("hsa-miR195", "hsa-miR296","hsa-miR497"))


###################################################
### code chunk number 19: functions
###################################################

# pulling out the target genes for each miRNAs
foo = function(y,x) match(y$transcript_cluster_id, rownames(x))


# organizing "full" dataset of gene expression data into the dataframes of patients for each genes, with
# labels for separate tissue
org = function(df1,df2){
	m = foo(df2,df1)
	df = df1[m,]
	s = split(df1, rownames(df1))
	s = s[rownames(df)]
	names(s) = df2$hgnc_symbol
	ind = which(metadata$TissueType == "normal_mucosa")
	d = rep(names(data),sapply(data,function(x) nrow(x)/2))
	d = d[(d=="nonpoly") | (d=="poly")]
	lapply(s, function(u) data.frame(label = d, nm = u[ind], l = u[-ind]))
	
}



###################################################
### code chunk number 20: preparing_for_plotting
###################################################
indices = mapply(foo, targets, dm.tissue)

exp = lapply(indices, function(u) {
	lapply(1:length(data),function(v) {
		   dm.tissue[[v]][u,]
		   })
})

genes = lapply(exp,"[[",1)
limma = storage[[1]]
final = mapply(org,genes,limma)

np.p = storage[-1]
pn = vector("list",3) ;names(pn) = names(f)
for(i in 1:length(limma)){
	g = lapply(lapply(np.p,"[[",i), function (x){
			   rownames(x) = x$hgnc_symbol
			   x[limma[[i]]$hgnc_symbol,]
			   })
	
		pn[[i]] = g
	
}



###################################################
### code chunk number 21: plotting
###################################################
pdf(file = "full_plot_miRNA.pdf",width=15,height=10)
par(mfrow = c(2,2))
for(i in 1:length(final)){
	for(j in 1:length(final[[i]])){
		matplot(c(1,2),rbind(unlist(final[[i]][[j]][2]), unlist(final[[i]][[j]][3])),
				type="l", lwd=3, xlim = c(1.0, 2.17), 
				col = final[[i]][[j]]$label,
				lty = 1,
				xlab = " ", ylab = " ",xaxt = "n")
		title(main = paste(limma[[i]]$hgnc_symbol[j], names(f)[i],sep = "/"),
			  sub = paste("adjusted p-values:",
						  "full =", sprintf("%3.2e", limma[[i]][j,10]),
						  "poly =", sprintf("%3.2e", pn[[i]][[2]][j,10]), 
						  "npoly =", sprintf("%3.2e", pn[[i]][[1]][j,10]), sep = " "))
		axis(1,1:2, c("normal mucosa", "lesion"))	
		legend("bottomright", names(np.p), lty=c(1,1), lwd=c(2.5,2.5), col = c("black","red"))
		
	}
}



###################################################
### code chunk number 22: dev_off3
###################################################
y <- dev.off()


###################################################
### code chunk number 23: enviro
###################################################
sessionInfo()


