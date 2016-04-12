#Geo ID - http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4761
#GDS4761
#Probe Ids
#100157844_TGI_at - CLDN1 
#100128434_TGI_at - CLDN25

require(GEOquery)
gds <- getGEO("GDS4761")
eset <- GDS2eSet(gds)

require(Biobase)
dat <- exprs(eset)
selectedProbeIDs[1] = "100157844_TGI_at"
selectedProbeIDs[2] = "100128434_TGI_at"
naProbeSets <- apply(dat, 1, function(x) {sum(is.na(x))})
sum(naProbeSets != 0)
datNoNa <- dat[naProbeSets == 0, ]
df <- as.data.frame(datNoNa)
diffExpP <- c()
diffExpP = apply(df,1,function(x){
    class1Exp <- x[pData(eset)$tissue == "liver metastasis"]
    class2Exp <- x[pData(eset)$tissue == "local metastasis in the breast"]
    t.test(class1Exp, class2Exp)$p.value
})

heatmap(as.matrix(df[diffExpP < 0.00001,]))