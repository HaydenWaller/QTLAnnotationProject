### parkoh ###

## parkoh has been coded in the wrong way, so L. kohalensis allele are B and L. paranigra alleles are A


parkoh<-read.cross(format="csv", file="ParKoh_Comprehensive_CrossFile.csv", estimate.map=FALSE)
parkoh_jitter<-jittermap(parkoh)
parkoh_jitter<-sim.geno(parkoh_jitter, step=1, n.draws=1000, error.prob=0.001)
parkoh_jitter<-calc.genoprob(parkoh_jitter, step=1, error.prob=0.001)

## scanone

qtl_parkoh_imp<-scanone(parkoh_jitter,method="imp", pheno.col="pr")
qtl_parkoh_hk<-scanone(parkoh_jitter,method="hk", pheno.col="pr")
qtl_parkoh_ehk<-scanone(parkoh_jitter,method="ehk", pheno.col="pr")
qtl_parkoh_em<-scanone(parkoh_jitter,method="em", pheno.col="pr")

qtl_parkoh_imp_1000perm<-scanone(parkoh_jitter,method="imp", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE) # this takes about 2 minutes
qtl_parkoh_hk_1000perm<-scanone(parkoh_jitter,method="hk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
qtl_parkoh_em_1000perm<-scanone(parkoh_jitter,method="em", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
#qtl_parkoh_ehk_1000perm<-scanone(parkoh_jitter,method="ehk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)


## cim

#cim_parkoh_imp<-cim(parkoh_jitter, pheno.col="pr", method="imp", n.marcovar=3, window=10)
cim_parkoh_hk<-cim(parkoh_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20)
#cim_parkoh_ehk<-cim(parkoh_jitter, pheno.col="pr", method="ehk", n.marcovar=3, window=10)
#cim_parkoh_em<-cim(parkoh_jitter, pheno.col="pr", method="em", n.marcovar=3, window=10)

cim_parkoh_hk_1000perm<-cim(parkoh_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20, n.perm=1000)

#cim_parkoh_imp_w5<-cim(parkoh_jitter, pheno.col="pr", method="imp", n.marcovar=3, window=5)
#cim_parkoh_hk_w5<-cim(parkoh_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=5)
#cim_parkoh_ehk_w5<-cim(parkoh_jitter, pheno.col="pr", method="ehk", n.marcovar=3, window=5)
#cim_parkoh_em_w5<-cim(parkoh_jitter, pheno.col="pr", method="em", n.marcovar=3, window=5)

cairo_pdf("cim_parkoh_hk_5_20.pdf", width=8, height=3.5)
plot(cim_parkoh_hk, ylim=c(0,50));abline(h=summary(cim_parkoh_hk_1000perm)[1,1])
dev.off()

qtl_parkoh<-c("1","2","3","4","5","6","X")

## scantwo
twoqtl_parkoh_imp<-scantwo(parkoh_jitter, pheno.col="pr", method="imp", chr=qtl_parkoh, clean.output=TRUE)
twoqtl_parkoh_hk<-scantwo(parkoh_jitter, pheno.col="pr", method="hk", chr=qtl_parkoh, clean.output=TRUE)
twoqtl_parkoh_em<-scantwo(parkoh_jitter, pheno.col="pr", method="em", chr=qtl_parkoh, clean.output=TRUE)

# check for potential interactions:

plot(twoqtl_parkoh_imp, lower="cond-int")
plot(twoqtl_parkoh_hk, lower="cond-int")
plot(twoqtl_parkoh_em, lower="cond-int")

# potential interaction at c(2,4), c(2,5), c(4,5)

# check for secondary peaks at specific chromosomes:

plot(twoqtl_parkoh_imp, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_imp, chr=2, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_imp, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_imp, chr=4, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_imp, chr=5, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_imp, chr=6, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_imp, chr="X", lower="cond-int", upper="cond-add")


plot(twoqtl_parkoh_hk, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_hk, chr=2, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_hk, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_hk, chr=4, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_hk, chr=5, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_hk, chr=6, lower="cond-int", upper="cond-add")
plot(twoqtl_parkoh_hk, chr="X", lower="cond-int", upper="cond-add")


# permutate the scantwo function 
twoqtl_parkoh_hk_1000perm<-scantwo(parkoh_jitter, pheno.col="pr", method="hk", chr=qtl_parkoh, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes
summary(twoqtl_parkoh_hk, threshold=as.matrix(as.data.frame(summary(twoqtl_parkoh_hk_10000perm, 0.05)[1:5])[1,]))
penalties=calc.penalties(twoqtl_parkoh_hk_1000perm, alpha=0.05)
    main    heavy    light 
3.415696 5.533016 3.307139

max(qtl_parkoh_imp)
               chr pos  lod
S002565_116928   5  31 17.6


parkoh_init_qtl<-makeqtl(parkoh_jitter,chr=5, pos=31)
summary(fitqtl(parkoh_jitter, qtl=parkoh_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1))
parkoh_init_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_init_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1)

parkoh_expand_qtl<-addtoqtl(parkoh_jitter,parkoh_init_qtl, 2, 96)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
parkoh_expand_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand_qtl)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
parkoh_expand_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand_qtl, method="imp", formula=  y ~ Q1 + Q2)
parkoh_expand_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2)

parkoh_expand2_qtl<-addtoqtl(parkoh_jitter,parkoh_expand_qtl, 4, 41.9)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
parkoh_expand2_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand2_qtl)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
parkoh_expand2_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand2_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3)
parkoh_expand2_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand2_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3)

parkoh_expand3_qtl<-addtoqtl(parkoh_jitter,parkoh_expand2_qtl, 1, 91.4)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4))
parkoh_expand3_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand3_qtl)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4))
parkoh_expand3_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand3_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4)
parkoh_expand3_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand3_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q2:Q4)) # pLOD is not sufficient to allow interaction to be added

parkoh_expand4_qtl<-addtoqtl(parkoh_jitter,parkoh_expand3_qtl, "X", 27)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))
parkoh_expand4_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand4_qtl)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))
parkoh_expand4_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand4_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5)
parkoh_expand4_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand4_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q2:Q4)) # this time it is about enough

parkoh_expand5_qtl<-addtoqtl(parkoh_jitter,parkoh_expand4_qtl, 2, 200)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q2:Q4))
parkoh_expand5_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand5_qtl)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)) # when you drop the interaction, it no longer satisfies the pLOD criteria, so drop it!
parkoh_expand5_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand5_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)
parkoh_expand5_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand5_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)

parkoh_expand6_qtl<-addtoqtl(parkoh_jitter,parkoh_expand5_qtl, 6, 57)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7))
parkoh_expand6_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand6_qtl)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7))
parkoh_expand6_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand6_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7)
parkoh_expand6_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand6_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7)

parkoh_expand7_qtl<-addtoqtl(parkoh_jitter,parkoh_expand6_qtl, 3, 126.8)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8))
parkoh_expand7_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand7_qtl)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8))
parkoh_expand7_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand7_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8)
parkoh_expand7_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand7_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8)


parkoh_expand8_qtl<-addtoqtl(parkoh_jitter,parkoh_expand7_qtl, 3, 43.0)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9))
parkoh_expand8_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand8_qtl)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9))
parkoh_expand8_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand8_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9)
parkoh_expand8_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand8_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9)
summary(fitqtl(parkoh_jitter, qtl=parkoh_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9))

cairo_pdf("ParKoh_fitqtl.pdf", width=8, height=3.5)
plotLodProfile(parkoh_expand8_qtl,showallchr=TRUE); abline(h=3.41)
dev.off()

## LOD 1 and LOD 2 intervals
rbind(
lodint(parkoh_expand8_qtl, qtl.index=1, drop=1, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=1, drop=2, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=2, drop=1, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=2, drop=2, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=3, drop=1, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=3, drop=2, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=4, drop=1, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=4, drop=2, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=5, drop=1, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=5, drop=2, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=6, drop=1, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=6, drop=2, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=7, drop=1, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=7, drop=2, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=8, drop=1, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=8, drop=2, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=9, drop=1, expandtomarker=TRUE),
lodint(parkoh_expand8_qtl, qtl.index=9, drop=2, expandtomarker=TRUE))


parkoh_QTLintervals<-data.frame(rbind(
cbind(paste("LODint1"),as.character(lodint(parkoh_expand8_qtl, qtl.index=1, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=1, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=1, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(parkoh_expand8_qtl, qtl.index=1, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=1, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=1, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(parkoh_expand8_qtl, qtl.index=2, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=2, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=2, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(parkoh_expand8_qtl, qtl.index=2, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=2, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=2, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(parkoh_expand8_qtl, qtl.index=3, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=3, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=3, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(parkoh_expand8_qtl, qtl.index=3, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=3, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=3, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(parkoh_expand8_qtl, qtl.index=4, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=4, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=4, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(parkoh_expand8_qtl, qtl.index=4, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=4, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=4, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(parkoh_expand8_qtl, qtl.index=5, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=5, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=5, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(parkoh_expand8_qtl, qtl.index=5, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=5, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=5, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(parkoh_expand8_qtl, qtl.index=6, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=6, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=6, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(parkoh_expand8_qtl, qtl.index=6, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=6, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=6, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(parkoh_expand8_qtl, qtl.index=7, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=7, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=7, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(parkoh_expand8_qtl, qtl.index=7, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=7, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=7, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(parkoh_expand8_qtl, qtl.index=8, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=8, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=8, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(parkoh_expand8_qtl, qtl.index=8, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=8, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=8, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(parkoh_expand8_qtl, qtl.index=9, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=9, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=9, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(parkoh_expand8_qtl, qtl.index=9, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=9, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=9, drop=2, expandtomarker=TRUE))))]))
)

parkoh_QTLintervals<-data.frame(rbind(
cbind(paste("LODint1.5"),as.character(lodint(parkoh_expand8_qtl, qtl.index=2, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=2, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=2, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkoh_expand8_qtl, qtl.index=6, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=6, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=6, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkoh_expand8_qtl, qtl.index=8, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=8, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=8, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkoh_expand8_qtl, qtl.index=9, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=9, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=9, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkoh_expand8_qtl, qtl.index=4, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=4, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=4, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkoh_expand8_qtl, qtl.index=3, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=3, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=3, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkoh_expand8_qtl, qtl.index=1, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=1, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=1, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkoh_expand8_qtl, qtl.index=7, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=7, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=7, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkoh_expand8_qtl, qtl.index=5, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkoh_expand8_qtl, qtl.index=5, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkoh_expand8_qtl, qtl.index=5, drop=1.5, expandtomarker=TRUE))))]))
)


parkoh_QTLintervals$pos<-0
for(i in 1:nrow(parkoh_QTLintervals)) { parkoh_QTLintervals$pos[i]=find.markerpos(parkoh_jitter,as.character(parkoh_QTLintervals$X3[i]))$pos}

# table with QTL effects
parkoh_QTLeffects<-data.frame(rbind(
data.frame(marker=find.marker(parkoh_jitter, 2,70.0),markerpos=find.markerpos(parkoh_jitter,as.character(find.marker(parkoh_jitter, 2,70.0)))$pos,matrix(effectplot(parkoh_jitter,"2@70.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkoh_jitter, 2,198.1),markerpos=find.markerpos(parkoh_jitter,as.character(find.marker(parkoh_jitter, 2,198.1)))$pos,matrix(effectplot(parkoh_jitter,"2@7198.1", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkoh_jitter, 3,74.0),markerpos=find.markerpos(parkoh_jitter,as.character(find.marker(parkoh_jitter, 3,74.0)))$pos,matrix(effectplot(parkoh_jitter,"3@74.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkoh_jitter, 3,165.6),markerpos=find.markerpos(parkoh_jitter,as.character(find.marker(parkoh_jitter, 3,165.6)))$pos,matrix(effectplot(parkoh_jitter,"3@165.6", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkoh_jitter, 1,86.3),markerpos=find.markerpos(parkoh_jitter,as.character(find.marker(parkoh_jitter, 1,86.3)))$pos,matrix(effectplot(parkoh_jitter,"1@86.3", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkoh_jitter, 4,39.0),markerpos=find.markerpos(parkoh_jitter,as.character(find.marker(parkoh_jitter, 4,39.0)))$pos,matrix(effectplot(parkoh_jitter,"4@39.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkoh_jitter, 5,31.0),markerpos=find.markerpos(parkoh_jitter,as.character(find.marker(parkoh_jitter, 5,31.0)))$pos,matrix(effectplot(parkoh_jitter,"5@31.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkoh_jitter, 6,57.0),markerpos=find.markerpos(parkoh_jitter,as.character(find.marker(parkoh_jitter, 6,57.0)))$pos,matrix(effectplot(parkoh_jitter,"6@57.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkoh_jitter, "X",41.0),markerpos=find.markerpos(parkoh_jitter,as.character(find.marker(parkoh_jitter, "X",41.0)))$pos,matrix(effectplot(parkoh_jitter,"X@41.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)))
)


parkoh_QTLeffects$est<-(parkoh_QTLeffects$X1-parkoh_QTLeffects$X3)/(-2)
parkoh_QTLeffects$est[8]<-(parkoh_QTLeffects$X1[9]-parkoh_QTLeffects$X2[9])/(-2)
parkoh_QTLeffects$PVE<-100*parkoh_QTLeffects$est/3.7
parkoh_QTLeffects$QTL<-c("2a","2b","3a","3b","1","4","5a","6","X")


cairo_pdf("ParKoh_effectplots.pdf", height=28, width=3)
par(mfrow = c(9,1))
effectplot(parkoh_jitter,"2@70.0", pheno.col="pr", ylim=c(1,3))
effectplot(parkoh_jitter,"2@198.1", pheno.col="pr", ylim=c(1,3))
effectplot(parkoh_jitter,"3@74.0", pheno.col="pr", ylim=c(1,3))
effectplot(parkoh_jitter,"3@165.6", pheno.col="pr", ylim=c(1,3))
effectplot(parkoh_jitter,"1@86.3", pheno.col="pr", ylim=c(1,3))
effectplot(parkoh_jitter,"4@39.0", pheno.col="pr", ylim=c(1,3))
effectplot(parkoh_jitter,"5@31.0", pheno.col="pr", ylim=c(1,3))
effectplot(parkoh_jitter,"6@57.0", pheno.col="pr", ylim=c(1,3))
effectplot(parkoh_jitter,"X@41.0", pheno.col="pr", ylim=c(1,3))
dev.off()


# Bayesian confidence intervals
parkoh_QTLbayesintervals<-data.frame(locus=character(), position=numeric())
for(i in 1:9) {
	temp_df<-data.frame(locus=names(pull.map(parkoh_jitter)[[as.character(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]),position=matrix(pull.map(parkoh_jitter)[[as.character(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]), LG=as.character(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])
) 

	from=min(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	to=max(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)

	parkoh_QTLbayesintervals<-rbind(parkoh_QTLbayesintervals,temp_df[max(which(temp_df$position<=from)):min(which(temp_df$position>=to)),])
	}
	
		
parkoh_QTLbayesintervals$scaffold<-gsub("_.*","",parkoh_QTLbayesintervals$locus)

write.table(parkoh_QTLbayesintervals,"parkoh_QTLbayesintervals.txt", quote=FALSE, row.names=FALSE, sep="\t")
