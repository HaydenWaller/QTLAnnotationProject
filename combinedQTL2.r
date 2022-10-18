## mean pulse rate Lkoh = 3.7, Lpar=0.7, Lkon = 1.8-2.3, Lpru 1.5-2.0
parkoh_effects<-data.frame(QTL=c(1,2,2,3,3,4,5,6,"X"),effectsize=c(0.1822,0.2338,0.09,0.0664,0.0505,0.1990,0.2679,0.0705,0.0797))
parkoh_effects$releffectsize=parkoh_effects$effectsize/(3.7-0.7)*100

library(qtl)

qtl_parkoh<-c("1","2","4","5","6","X")
twoqtl_parkoh_hk_1000perm<-scantwo(parkoh_jitter, pheno.col="pr", method="hk", chr=qtl_parkoh, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes

qtl_parkon<-c("1","2","3","4","5","6","X")
twoqtl_parkon_hk_1000perm<-scantwo(parkon_jitter, pheno.col="pr", method="hk", chr=qtl_parkon, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 1,5 hours

qtl_prukoh<-c("1","2","3","4","5","6","X")
twoqtl_prukoh_hk_1000perm<-scantwo(prukoh_jitter, pheno.col="pr", method="hk", chr=qtl_prukoh, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes


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

cairo_pdf("GRAD_SKEWL/Projects_Experiments/RNA-seq/reanalysis_project/cim_parkoh_hk_5_20.pdf", width=8, height=3.5)
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

cairo_pdf("GRAD_SKEWL/Projects_Experiments/RNA-seq/reanalysis_project/figures/ParKoh_fitqtl.pdf", width=8, height=3.5)
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


cairo_pdf("/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/RNA-seq/reanalysis_project/ParKoh_effectplots.pdf", height=28, width=3)
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
CIsize_parkoh<-c()
for(i in 1:9) {
	temp_df<-data.frame(locus=names(pull.map(parkoh_jitter)[[as.character(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]),position=matrix(pull.map(parkoh_jitter)[[as.character(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]), LG=as.character(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])
) 

	from=min(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	to=max(bayesint(parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	CIsize_parkoh<-c(CIsize_parkoh,to-from)

	parkoh_QTLbayesintervals<-rbind(parkoh_QTLbayesintervals,temp_df[max(which(temp_df$position<=from)):min(which(temp_df$position>=to)),])
	}
	
		
parkoh_QTLbayesintervals$scaffold<-gsub("_.*","",parkoh_QTLbayesintervals$locus)

write.table(parkoh_QTLbayesintervals,"/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/RNA-seq/reanalysis_project/parkoh_QTLbayesintervals.txt", quote=FALSE, row.names=FALSE, sep="\t")

additive_effects_parkoh<-matrix(summary(fitqtl(parkoh_jitter, qtl=parkoh_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9))$ests[c(2,4,6,8,10,11,13,15,17),1])
summary(lm(CIsize_parkoh~additive_effects_parkoh))


parkon_QTLintervals$pos<-0
for(i in 1:nrow(parkon_QTLintervals)) { parkon_QTLintervals$pos[i]=find.markerpos(parkon_jitter,as.character(parkon_QTLintervals$X3[i]))$pos}

peakmarkerlist_parkon=c(substr(as.matrix(find.flanking(parkon_jitter,1,79.1)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,2,12.7)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,2,58.1)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,3,60.0)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,4,48.0)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,5,32.5)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,5,57.6)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,6,47.0)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,"X",55.0)),0,7))

write.table(peakmarkerlist_parkon,"/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/RNA-seq/reanalysis_project/peakmarkerlist_parkon.txt", quote=FALSE, row.names=FALSE, sep="\t")

fitdistr(parkoh_QTLeffects$est, "exponential") # fit exponential distribution to data
cairo_pdf("/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/RNA-seq/reanalysis_project/parkoh_effectsize_distr.pdf", width=3, height=3)
plot(density(rexp(10000, rate=fitdistr(parkoh_QTLeffects$est, "exponential")$estimate)),xlim=c(0,1))
hist(parkoh_QTLeffects$est, prob=TRUE, breaks=4, add=TRUE)
dev.off()

parkoh_QTL2LODintervals<-data.frame(locus=character(), position=numeric())
for(i in 1:9) {
	temp_df<-data.frame(locus=names(pull.map(parkoh_jitter)[[as.character(lodint(drop=2.0,parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]),position=matrix(pull.map(parkoh_jitter)[[as.character(lodint(drop=2.0,parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]), LG=as.character(lodint(drop=2.0,parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])
) 

	from=min(lodint(drop=2.0,parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	to=max(lodint(drop=2.0,parkoh_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)

	parkoh_QTL2LODintervals<-rbind(parkoh_QTL2LODintervals,temp_df[max(which(temp_df$position<=from)):min(which(temp_df$position>=to)),])
	}
	
	#candidate genes location on the map
	parkoh_no_marker_names<-parkoh
names(parkoh_no_marker_names$geno$"1"$map)<-paste0("1_",as.numeric(as.factor(names(parkoh_no_marker_names$geno$"1"$map))))
names(parkoh_no_marker_names$geno$"2"$map)<-paste0("2_",as.numeric(as.factor(names(parkoh_no_marker_names$geno$"2"$map))))
names(parkoh_no_marker_names$geno$"3"$map)<-paste0("3_",as.numeric(as.factor(names(parkoh_no_marker_names$geno$"3"$map))))
names(parkoh_no_marker_names$geno$"4"$map)<-paste0("4_",as.numeric(as.factor(names(parkoh_no_marker_names$geno$"4"$map))))
names(parkoh_no_marker_names$geno$"5"$map)<-paste0("5_",as.numeric(as.factor(names(parkoh_no_marker_names$geno$"5"$map))))
names(parkoh_no_marker_names$geno$"6"$map)<-paste0("6_",as.numeric(as.factor(names(parkoh_no_marker_names$geno$"6"$map))))
names(parkoh_no_marker_names$geno$"7"$map)<-paste0("7_",as.numeric(as.factor(names(parkoh_no_marker_names$geno$"7"$map))))
names(parkoh_no_marker_names$geno$"X"$map)<-paste0("X_",as.numeric(as.factor(names(parkoh_no_marker_names$geno$"X"$map))))

colnames(parkoh_no_marker_names$geno$"1"$data)<-paste0("1_",as.numeric(as.factor(colnames(parkoh_no_marker_names$geno$"1"$data))))
colnames(parkoh_no_marker_names$geno$"2"$data)<-paste0("2_",as.numeric(as.factor(colnames(parkoh_no_marker_names$geno$"2"$data))))
colnames(parkoh_no_marker_names$geno$"3"$data)<-paste0("3_",as.numeric(as.factor(colnames(parkoh_no_marker_names$geno$"3"$data))))
colnames(parkoh_no_marker_names$geno$"4"$data)<-paste0("4_",as.numeric(as.factor(colnames(parkoh_no_marker_names$geno$"4"$data))))
colnames(parkoh_no_marker_names$geno$"5"$data)<-paste0("5_",as.numeric(as.factor(colnames(parkoh_no_marker_names$geno$"5"$data))))
colnames(parkoh_no_marker_names$geno$"6"$data)<-paste0("6_",as.numeric(as.factor(colnames(parkoh_no_marker_names$geno$"6"$data))))
colnames(parkoh_no_marker_names$geno$"7"$data)<-paste0("7_",as.numeric(as.factor(colnames(parkoh_no_marker_names$geno$"7"$data))))
colnames(parkoh_no_marker_names$geno$"X"$data)<-paste0("X_",as.numeric(as.factor(colnames(parkoh_no_marker_names$geno$"X"$data))))

names(parkoh_no_marker_names$geno$"1"$map)[65]<-"candidate_gene___ato"
colnames(parkoh_no_marker_names$geno$"1"$data)[65]<-"candidate_gene___ato"

names(parkoh_no_marker_names$geno$"1"$map)[67]<-"candidate_gene___e"
colnames(parkoh_no_marker_names$geno$"1"$data)[67]<-"candidate_gene___e"

names(parkoh_no_marker_names$geno$"1"$map)[59]<-"candidate_gene___slo"
colnames(parkoh_no_marker_names$geno$"1"$data)[59]<-"candidate_gene___slo"


names(parkoh_no_marker_names$geno$"3"$map)[37]<-"candidate_gene___btv1"
colnames(parkoh_no_marker_names$geno$"3"$data)[37]<-"candidate_gene___btv1"

names(parkoh_no_marker_names$geno$"3"$map)[40]<-"candidate_gene___btv2"
colnames(parkoh_no_marker_names$geno$"3"$data)[40]<-"candidate_gene___btv2"

names(parkoh_no_marker_names$geno$"3"$map)[29]<-"candidate_gene___fru"
colnames(parkoh_no_marker_names$geno$"3"$data)[29]<-"candidate_gene___fru"

names(parkoh_no_marker_names$geno$"3"$map)[34]<-"candidate_gene___fru2"
colnames(parkoh_no_marker_names$geno$"3"$data)[34]<-"candidate_gene___fru2"


names(parkoh_no_marker_names$geno$"4"$map)[7]<-"candidate_gene___cac"
colnames(parkoh_no_marker_names$geno$"4"$data)[7]<-"candidate_gene___cac"

names(parkoh_no_marker_names$geno$"4"$map)[8]<-"candidate_gene___cac2"
colnames(parkoh_no_marker_names$geno$"4"$data)[8]<-"candidate_gene___cac2"


names(parkoh_no_marker_names$geno$"X"$map)[4]<-"candidate_gene___dxs"
colnames(parkoh_no_marker_names$geno$"X"$data)[4]<-"candidate_gene___dxs"

names(parkoh_no_marker_names$geno$"X"$map)[5]<-"candidate_gene___dxs2"
colnames(parkoh_no_marker_names$geno$"X"$data)[5]<-"candidate_gene___dxs2"


names(parkoh_no_marker_names$geno$"X"$map)[21]<-"candidate_gene___syn_1"
colnames(parkoh_no_marker_names$geno$"X"$data)[21]<-"candidate_gene___syn_1"

names(parkoh_no_marker_names$geno$"X"$map)[22]<-"candidate_gene___syn1_2"
colnames(parkoh_no_marker_names$geno$"X"$data)[22]<-"candidate_gene___syn1_2"



names(parkoh_no_marker_names$geno$"2"$map)[81]<-"candidate_gene___per"
colnames(parkoh_no_marker_names$geno$"2"$data)[81]<-"candidate_gene___per"


parkoh_no_marker_names_qtl<-scanone(jittermap(parkoh_no_marker_names),pheno.col="pr",method="hk")

cairo_pdf("ParKoh_qtl_markernames.pdf", width=8, height=3.5)
plot(parkoh_no_marker_names_qtl, show.marker.names=TRUE)
dev.off()

#parkoh_init_qtl<-makeqtl(parkoh_jitter,chr=c(1,2,4,5,"X"), pos=c(113.2,69.2,35.8,31.0,36.8))
#summary(fitqtl(parkoh_jitter, qtl=parkoh_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))
#parkoh_init_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_init_qtl)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))

#parkoh_expand_qtl<-addtoqtl(parkoh_jitter,parkoh_init_qtl, 1, 47)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
#parkoh_expand_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand_qtl)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)) # not significant


#parkoh_expand_qtl<-addtoqtl(parkoh_jitter,parkoh_init_qtl, 2, 113)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
#parkoh_expand_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand_qtl)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))


#parkoh_expand_qtl<-addtoqtl(parkoh_jitter,parkoh_init_qtl, 4, 70)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
#parkoh_expand_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand_qtl)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)) NOT significant


#parkoh_expand_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q1:Q2))
#parkoh_expand_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)

#parkoh_expand2_qtl<-addtoqtl(parkoh_jitter,parkoh_expand_qtl, c(3,6), c(137, 57))
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q2))
#parkoh_expand2_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand2_qtl)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q2))

#parkoh_expand2_qtl_addint<-addint(parkoh_jitter, pheno.col="pr",  qtl=parkoh_expand2_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q2)
#parkoh_expand2_qtl_add <- addqtl(parkoh_jitter, qtl=parkoh_expand2_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q2)

#parkoh_expand3_qtl<-addtoqtl(parkoh_jitter,parkoh_expand2_qtl, 3, 23.9)
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2))
#parkoh_expand3_qtl<-refineqtl(parkoh_jitter, pheno.col="pr", method="imp", qtl = parkoh_expand3_qtl, )
#summary(fitqtl(parkoh_jitter, qtl=parkoh_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2))

#rbind(
#bayesint(parkoh_expand3_qtl, qtl.index=1),
#bayesint(parkoh_expand3_qtl, qtl.index=2),
#bayesint(parkoh_expand3_qtl, qtl.index=3),
####bayesint(parkoh_expand3_qtl, qtl.index=4),
#bayesint(parkoh_expand3_qtl, qtl.index=5),
#bayesint(parkoh_expand3_qtl, qtl.index=6),
#bayesint(parkoh_expand3_qtl, qtl.index=7),
###bayesint(parkoh_expand3_qtl, qtl.index=8),
#bayesint(parkoh_expand3_qtl, qtl.index=9))

## multiple QTL models / stepwiseqtl
#summary(fitqtl(parkoh_jitter_sim, qtl=parkoh_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr"))
#parkoh_add_qtl <- addqtl(parkoh_jitter_sim, qtl=parkoh_init_qtl, pheno.col="pr", method="imp") # check for additional QTL after correcting for main QTL


#stepwise_qtl_parkoh_imp<-stepwiseqtl(parkoh_jitter_sim, method="imp", model="normal", scan.pairs=FALSE, pheno.col="pr", penalties=penalties, keeptrace=TRUE)
#stepwise_qtl_parkoh_hk<-stepwiseqtl(parkoh_jitter_prob, method="hk", model="normal", scan.pairs=FALSE, pheno.col="pr", penalties=penalties, keeptrace=TRUE)

#stepwise_qtl_parkoh_imp_dropQ9<-dropfromqtl(stepwise_qtl_parkoh_imp, 9) # Q9 is unlikely, drop from the model
#stepwise_qtl_parkoh_imp_dropQ9<-refineqtl(parkoh_jitter_sim,qtl=stepwise_qtl_parkoh_imp_dropQ9, method="imp", pheno.col="pr") # refine locations

#stepwise_qtl_parkoh_imp_dropQ9_addint<-addint(parkoh_jitter_sim, pheno.col="pr", qtl=stepwise_qtl_parkoh_imp_dropQ9, method="imp") # check for interactions
#stepwise_qtl_parkoh_imp_final<-refineqtl(parkoh_jitter_sim, pheno.col="pr", method="imp", qtl = stepwise_qtl_parkoh_imp_dropQ9, formula = y~Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q1:Q2+Q1:Q7) # refine


### parkon ####

# parkon has been coded the wrong way, L. paranigra alleles are B, L. kona alles are A

parkon<-read.cross(format="csv", file="ParKon_Comprehensive_CrossFile.csv", estimate.map=FALSE)
parkon_jitter<-jittermap(parkon)
parkon_jitter<-sim.geno(parkon_jitter, step=1, n.draws=1000, error.prob=0.001)
parkon_jitter<-calc.genoprob(parkon_jitter, step=1, error.prob=0.001)

## scanone

qtl_parkon_imp<-scanone(parkon_jitter,method="imp", pheno.col="pr")
qtl_parkon_hk<-scanone(parkon_jitter,method="hk", pheno.col="pr")
qtl_parkon_ehk<-scanone(parkon_jitter,method="ehk", pheno.col="pr")
qtl_parkon_em<-scanone(parkon_jitter,method="em", pheno.col="pr")

qtl_parkon_imp_1000perm<-scanone(parkon_jitter,method="imp", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE) # this takes about 2 minutes
qtl_parkon_hk_1000perm<-scanone(parkon_jitter,method="hk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
qtl_parkon_em_1000perm<-scanone(parkon_jitter,method="em", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
#qtl_parkon_ehk_1000perm<-scanone(parkon_jitter,method="ehk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)


## cim

#cim_parkon_imp<-cim(parkon_jitter, pheno.col="pr", method="imp", n.marcovar=3, window=10)
cim_parkon_hk<-cim(parkon_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20)
#cim_parkon_ehk<-cim(parkon_jitter, pheno.col="pr", method="ehk", n.marcovar=3, window=10)
#cim_parkon_em<-cim(parkon_jitter, pheno.col="pr", method="em", n.marcovar=3, window=10)

cim_parkon_hk_1000perm<-cim(parkon_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20, n.perm=1000)

cairo_pdf("cim_parkon_hk_5_20.pdf", width=8, height=3.5)
plot(cim_parkon_hk, ylim=c(0,90));abline(h=summary(cim_parkon_hk_1000perm)[1,1])
dev.off()

qtl_parkon<-c("1","2","3","4","5","6","X")

## scantwo
twoqtl_parkon_imp<-scantwo(parkon_jitter, pheno.col="pr", method="imp", chr=qtl_parkon, clean.output=TRUE)
twoqtl_parkon_hk<-scantwo(parkon_jitter, pheno.col="pr", method="hk", chr=qtl_parkon, clean.output=TRUE)

plot(twoqtl_parkon_imp, lower="cond-int", upper="cond-add", c("3","5"))
plot(twoqtl_parkon_hk, lower="cond-int")

plot(twoqtl_parkon_imp, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_parkon_imp, chr=2, lower="cond-int", upper="cond-add")
# potential second qtl on chr 2
plot(twoqtl_parkon_imp, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_parkon_imp, chr=4, lower="cond-int", upper="cond-add")
plot(twoqtl_parkon_imp, chr=5, lower="cond-int", upper="cond-add")
# potential second qtl on chr 5
plot(twoqtl_parkon_imp, chr=6, lower="cond-int", upper="cond-add")
plot(twoqtl_parkon_imp, chr="X", lower="cond-int", upper="cond-add")

# permutate the scantwo function for LOD  penalties
twoqtl_parkon_hk_1000perm<-scantwo(parkon_jitter, pheno.col="pr", method="hk", chr=qtl_parkon, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 1,5 hours
summary(twoqtl_parkon_hk, threshold=as.matrix(as.data.frame(summary(twoqtl_parkon_hk_10000perm, 0.05)[1:5])[1,]))
penalties=calc.penalties(twoqtl_parkon_hk_1000perm, alpha=0.05)

    main    heavy    light 
3.464517 5.361266 3.197938 

max(qtl_parkon_hk)
               chr  pos  lod
S004383_170239   5 32.5 40.3

parkon_init_qtl<-makeqtl(parkon_jitter,chr=5, pos=32.5)
summary(fitqtl(parkon_jitter, qtl=parkon_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1))
parkon_init_qtl_add <- addqtl(parkon_jitter, qtl=parkon_init_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1)

parkon_expand_qtl<-addtoqtl(parkon_jitter,parkon_init_qtl, 4, 48)
summary(fitqtl(parkon_jitter, qtl=parkon_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
parkon_expand_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand_qtl)
summary(fitqtl(parkon_jitter, qtl=parkon_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
parkon_expand_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand_qtl, method="imp", formula=  y ~ Q1 + Q2)
parkon_expand_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2)


parkon_expand2_qtl<-addtoqtl(parkon_jitter,parkon_expand_qtl, 2, 57.6)
summary(fitqtl(parkon_jitter, qtl=parkon_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
parkon_expand2_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand2_qtl)
summary(fitqtl(parkon_jitter, qtl=parkon_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
parkon_expand2_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand2_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3) # Q1:Q3 and Q2:Q3 are both potentially significant, forward and backward selection shows that adding both violates light penalty, Q1:Q3 increases the LOD more so keep that one
parkon_expand2_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand2_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q1:Q3)

parkon_expand3_qtl<-addtoqtl(parkon_jitter,parkon_expand2_qtl, 1, 79.1)
summary(fitqtl(parkon_jitter, qtl=parkon_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q1:Q3))
parkon_expand3_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand3_qtl)
summary(fitqtl(parkon_jitter, qtl=parkon_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 +  Q1:Q3))
parkon_expand3_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand3_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 +  Q1:Q3) # now Q2:Q3 also is significant
parkon_expand3_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand3_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 +  Q1:Q3 + Q2:Q3)


parkon_expand4_qtl<-addtoqtl(parkon_jitter,parkon_expand3_qtl, 3, 60.7)
summary(fitqtl(parkon_jitter, qtl=parkon_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3))
parkon_expand4_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand4_qtl)
summary(fitqtl(parkon_jitter, qtl=parkon_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3))
parkon_expand4_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand4_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3)
parkon_expand4_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand4_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3 + Q1:Q5)
summary(fitqtl(parkon_jitter, qtl=parkon_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3 + Q1:Q5))

parkon_expand5_qtl<-addtoqtl(parkon_jitter,parkon_expand4_qtl, "X", 60)
summary(fitqtl(parkon_jitter, qtl=parkon_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q1:Q3 + Q2:Q3 + Q1:Q5))
parkon_expand5_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand5_qtl)
summary(fitqtl(parkon_jitter, qtl=parkon_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q1:Q3 + Q2:Q3 + Q1:Q5))
parkon_expand5_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand5_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q1:Q3 + Q2:Q3 + Q1:Q5)
parkon_expand5_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand5_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q1:Q3 + Q2:Q3 + Q1:Q5)


parkon_expand5b_qtl<-addtoqtl(parkon_jitter,parkon_expand5_qtl, 4, 50)
summary(fitqtl(parkon_jitter, qtl=parkon_expand5b_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5))
parkon_expand5b_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand5b_qtl)
summary(fitqtl(parkon_jitter, qtl=parkon_expand5b_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5))
parkon_expand5b_qtl_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand5b_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5)
parkon_expand5b_qtl_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand5b_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 	+ Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5)


parkon_expand6_qtl<-addtoqtl(parkon_jitter,parkon_expand5_qtl, 5, 57.6)
summary(fitqtl(parkon_jitter, qtl=parkon_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5))
parkon_expand6_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand6_qtl)
summary(fitqtl(parkon_jitter, qtl=parkon_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5))
parkon_expand6_qtl_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand6_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5)
parkon_expand6_qtl_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand6_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 	+ Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5)

parkon_expand7_qtl<-addtoqtl(parkon_jitter,parkon_expand6_qtl, 2, 12.7)
summary(fitqtl(parkon_jitter, qtl=parkon_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q3 + Q2:Q3 + Q1:Q5))
parkon_expand7_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand7_qtl)
summary(fitqtl(parkon_jitter, qtl=parkon_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q3 + Q2:Q3 + Q1:Q5)) # only Q2:Q3 is significant
summary(fitqtl(parkon_jitter, qtl=parkon_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q2 + Q2:Q3)) # Q1:Q2 is also significant
parkon_expand7_qtl_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand7_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q2:Q3) 
parkon_expand7_qtl_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand7_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 	+ Q6 + Q7 + Q8 + Q2:Q3 + Q1:Q2)

cairo_pdf("ParKon_fitqtl.pdf", width=8, height=3.5)
plotLodProfile(parkon_expand7_qtl,showallchr=TRUE, qtl.labels=FALSE); abline(h=3.46)
dev.off()

parkon_expand8_qtl<-addtoqtl(parkon_jitter,parkon_expand7_qtl, 6, 58.0)
summary(fitqtl(parkon_jitter, qtl=parkon_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3))
# strangely enough, refine here messes things up, lower LOD score, effect size, etc... So only refine QTL 1, which changes nothing (but is necessary for downstream stuff)
parkon_expand8_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand8_qtl, formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3)
summary(fitqtl(parkon_jitter, qtl=parkon_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3)) 
parkon_expand8_qtl_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand8_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3)  # Q3:Q7 is also significant
#parkon_expand8_qtl_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand8_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 	+ Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3)
summary(fitqtl(parkon_jitter, qtl=parkon_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3 + Q3:Q7))

cairo_pdf("ParKon_fitqtl_new.pdf", width=8, height=3.5)
plotLodProfile(parkon_expand8_qtl,showallchr=TRUE, qtl.labels=FALSE); abline(h=3.46)
dev.off()
#parkon_expand7_qtl<-addtoqtl(parkon_jitter,parkon_expand6_qtl, 5, 57.6)
#summary(fitqtl(parkon_jitter, qtl=parkon_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q3 + Q2:Q3 + Q1:Q5+ Q2:Q7))
#parkon_expand7_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand7_qtl)
#summary(fitqtl(parkon_jitter, qtl=parkon_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q3 + Q2:Q3 + Q1:Q5+ Q2:Q7))
#parkon_expand7_qtl_qtl_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand7_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q3 + Q2:Q3 + Q1:Q5+ Q2:Q7)
#parkon_expand7_qtl_qtl_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand7_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 	+ Q6 + Q7 + Q8 + Q1:Q3 + Q2:Q3 + Q1:Q5+ Q2:Q7) # all the ones added after this one are less than 1cM apart from existing peaks
#summary(fitqtl(parkon_jitter, qtl=parkon_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q3 + Q2:Q3 + Q1:Q5+ Q2:Q7 +Q1:Q8))

# parkon_QTLintervals<-data.frame(rbind(
# cbind(paste("LODint1"),as.character(lodint(parkon_expand7_qtl, qtl.index=1, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=1, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=1, drop=1, expandtomarker=TRUE))))]),
# cbind(paste("LODint2"),as.character(lodint(parkon_expand7_qtl, qtl.index=1, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=1, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=1, drop=2, expandtomarker=TRUE))))]),
# cbind(paste("LODint1"),as.character(lodint(parkon_expand7_qtl, qtl.index=2, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=2, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=2, drop=1, expandtomarker=TRUE))))]),
# cbind(paste("LODint2"),as.character(lodint(parkon_expand7_qtl, qtl.index=2, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=2, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=2, drop=2, expandtomarker=TRUE))))]),
# cbind(paste("LODint1"),as.character(lodint(parkon_expand7_qtl, qtl.index=3, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=3, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=3, drop=1, expandtomarker=TRUE))))]),
# cbind(paste("LODint2"),as.character(lodint(parkon_expand7_qtl, qtl.index=3, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=3, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=3, drop=2, expandtomarker=TRUE))))]),
# cbind(paste("LODint1"),as.character(lodint(parkon_expand7_qtl, qtl.index=4, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=4, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=4, drop=1, expandtomarker=TRUE))))]),
# cbind(paste("LODint2"),as.character(lodint(parkon_expand7_qtl, qtl.index=4, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=4, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=4, drop=2, expandtomarker=TRUE))))]),
# cbind(paste("LODint1"),as.character(lodint(parkon_expand7_qtl, qtl.index=5, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=5, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=5, drop=1, expandtomarker=TRUE))))]),
# cbind(paste("LODint2"),as.character(lodint(parkon_expand7_qtl, qtl.index=5, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=5, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=5, drop=2, expandtomarker=TRUE))))]),
# cbind(paste("LODint1"),as.character(lodint(parkon_expand7_qtl, qtl.index=6, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=6, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=6, drop=1, expandtomarker=TRUE))))]),
# cbind(paste("LODint2"),as.character(lodint(parkon_expand7_qtl, qtl.index=6, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=6, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=6, drop=2, expandtomarker=TRUE))))]),
# cbind(paste("LODint1"),as.character(lodint(parkon_expand7_qtl, qtl.index=7, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=7, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=7, drop=1, expandtomarker=TRUE))))]),
# cbind(paste("LODint2"),as.character(lodint(parkon_expand7_qtl, qtl.index=7, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=7, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=7, drop=2, expandtomarker=TRUE))))]),
# cbind(paste("LODint1"),as.character(lodint(parkon_expand7_qtl, qtl.index=8, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=8, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=8, drop=1, expandtomarker=TRUE))))]),
# cbind(paste("LODint2"),as.character(lodint(parkon_expand7_qtl, qtl.index=8, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand7_qtl, qtl.index=8, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand7_qtl, qtl.index=8, drop=2, expandtomarker=TRUE))))]))
# )


parkon_QTLintervals<-data.frame(rbind(
cbind(paste("LODint1.5"),as.character(lodint(parkon_expand8_qtl, qtl.index=8, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand8_qtl, qtl.index=8, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand8_qtl, qtl.index=8, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkon_expand8_qtl, qtl.index=2, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand8_qtl, qtl.index=2, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand8_qtl, qtl.index=2, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkon_expand8_qtl, qtl.index=5, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand8_qtl, qtl.index=5, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand8_qtl, qtl.index=5, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkon_expand8_qtl, qtl.index=4, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand8_qtl, qtl.index=4, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand8_qtl, qtl.index=4, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkon_expand8_qtl, qtl.index=1, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand8_qtl, qtl.index=1, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand8_qtl, qtl.index=1, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkon_expand8_qtl, qtl.index=7, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand8_qtl, qtl.index=7, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand8_qtl, qtl.index=7, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkon_expand8_qtl, qtl.index=6, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand8_qtl, qtl.index=6, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand8_qtl, qtl.index=6, drop=1.5, expandtomarker=TRUE))))]),
cbind(paste("LODint1.5"),as.character(lodint(parkon_expand8_qtl, qtl.index=9, drop=1.5, expandtomarker=TRUE)[1,1]),rownames(lodint(parkon_expand8_qtl, qtl.index=9, drop=1.5, expandtomarker=TRUE))[c(1,length(rownames(lodint(parkon_expand8_qtl, qtl.index=9, drop=1.5, expandtomarker=TRUE))))]))
)

parkon_QTLbayesintervals<-data.frame(LG=character(), locus=character(), position=numeric())
CIsize_parkon<-c()
for(i in 1:9) {
	temp_df<-data.frame(LG= bayesint(parkon_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1], locus=names(pull.map(parkon_jitter)[[as.character(bayesint(parkon_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]),position=matrix(pull.map(parkon_jitter)[[as.character(bayesint(parkon_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]))

	from=min(bayesint(parkon_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	to=max(bayesint(parkon_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	CIsize_parkon<-c(CIsize_parkon,to-from)

	parkon_QTLbayesintervals<-rbind(parkon_QTLbayesintervals,temp_df[max(which(temp_df$position<=from)):min(which(temp_df$position>=to)),])
	}
	
parkon_QTLbayesintervals$scaffold<-gsub("_.*","",parkon_QTLbayesintervals$locus)

write.table(parkon_QTLbayesintervals,"parkon_QTLbayesintervals.txt", quote=FALSE, row.names=FALSE, sep="\t")

additive_effects_parkon<-matrix(summary(fitqtl(parkon_jitter, qtl=parkon_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3 + Q3:Q7))$ests[c(2,4,6,8,10,12,13,15,17),1])
additive_effects_parkon<-(-1)*additive_effects_parkon
summary(lm(CIsize_parkon~additive_effects_parkon))

parkon_QTLintervals$pos<-0
for(i in 1:nrow(parkon_QTLintervals)) { parkon_QTLintervals$pos[i]=find.markerpos(parkon_jitter,as.character(parkon_QTLintervals$X3[i]))$pos}

peakmarkerlist_parkon=c(substr(as.matrix(find.flanking(parkon_jitter,1,79.1)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,2,12.7)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,2,58.1)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,3,60.0)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,4,48.0)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,5,32.5)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,5,57.6)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,6,47.0)),0,7),
substr(as.matrix(find.flanking(parkon_jitter,"X",55.0)),0,7))

write.table(peakmarkerlist_parkon,"peakmarkerlist_parkon.txt", quote=FALSE, row.names=FALSE, sep="\t")


parkon_QTLeffects<-data.frame(rbind(
data.frame(marker=find.marker(parkon_jitter, 2,58.1),markerpos=find.markerpos(parkon_jitter,as.character(find.marker(parkon_jitter, 2,58.1)))$pos,matrix(effectplot(parkon_jitter,"2@58.1", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkon_jitter, 2,12.7),markerpos=find.markerpos(parkon_jitter,as.character(find.marker(parkon_jitter, 2,12.7)))$pos,matrix(effectplot(parkon_jitter,"2@12.7", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkon_jitter, 3,60.0),markerpos=find.markerpos(parkon_jitter,as.character(find.marker(parkon_jitter, 3,60.0)))$pos,matrix(effectplot(parkon_jitter,"3@60.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkon_jitter, 1,79.1),markerpos=find.markerpos(parkon_jitter,as.character(find.marker(parkon_jitter, 1,79.1)))$pos,matrix(effectplot(parkon_jitter,"1@79.1", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkon_jitter, 4,48.0),markerpos=find.markerpos(parkon_jitter,as.character(find.marker(parkon_jitter, 4,48.0)))$pos,matrix(effectplot(parkon_jitter,"4@48.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkon_jitter, 5,32.5),markerpos=find.markerpos(parkon_jitter,as.character(find.marker(parkon_jitter, 5,32.5)))$pos,matrix(effectplot(parkon_jitter,"5@32.5", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkon_jitter, 5,57.6),markerpos=find.markerpos(parkon_jitter,as.character(find.marker(parkon_jitter, 5,57.6)))$pos,matrix(effectplot(parkon_jitter,"5@57.6", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkon_jitter, 6,47.0),markerpos=find.markerpos(parkon_jitter,as.character(find.marker(parkon_jitter, 5,47.0)))$pos,matrix(effectplot(parkon_jitter,"6@47.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(parkon_jitter, "X",55.0),markerpos=find.markerpos(parkon_jitter,as.character(find.marker(parkon_jitter, "X",55.0)))$pos,matrix(effectplot(parkon_jitter,"X@55.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)))
)

parkon_QTLeffects$est_single<-(parkon_QTLeffects$X1-parkon_QTLeffects$X3)/2
parkon_QTLeffects$est_single[9]<-(parkon_QTLeffects$X1[9]-parkon_QTLeffects$X2[9])/2
parkon_QTLeffects$est<-matrix(summary(fitqtl(parkon_jitter, qtl=parkon_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3))$est[c(4,15,10,8,6,2,13,17,12),1])[,1]*(-1)
parkon_QTLeffects$PVE<-100*parkon_QTLeffects$est/1.4
parkon_QTLeffects$QTL<-c("2a","2c","3a","1","4","5a","5b","6","X")

cairo_pdf("parkon_effectplots.pdf", height=28, width=3)
par(mfrow = c(9,1))
effectplot(parkon_jitter,"2@58.1", pheno.col="pr", ylim=c(1,3))
effectplot(parkon_jitter,"2@12.7", pheno.col="pr", ylim=c(1,3))
effectplot(parkon_jitter,"3@60.0", pheno.col="pr", ylim=c(1,3))
effectplot(parkon_jitter,"1@79.1", pheno.col="pr", ylim=c(1,3))
effectplot(parkon_jitter,"4@48.0", pheno.col="pr", ylim=c(1,3))
effectplot(parkon_jitter,"5@32.5", pheno.col="pr", ylim=c(1,3))
effectplot(parkon_jitter,"5@57.6", pheno.col="pr", ylim=c(1,3))
effectplot(parkon_jitter,"6@47.0", pheno.col="pr", ylim=c(1,3))
effectplot(parkon_jitter,"X@55.0", pheno.col="pr", ylim=c(1,3))
dev.off()


fitdistr(parkon_QTLeffects$est, "exponential") # fit exponential distribution to data
cairo_pdf("parkon_effectsize_distr.pdf", width=3, height=3)
plot(density(rexp(10000, rate=fitdistr(parkon_QTLeffects$est, "exponential")$estimate)),xlim=c(0,1))
hist(parkon_QTLeffects$est, prob=TRUE, add=TRUE)
dev.off()

#parkon_init_qtl<-makeqtl(parkon_jitter,chr=c(1,2,4,5,"X"), pos=c(67.0,64.2,48.0,32.5,66.7))
#summary(fitqtl(parkon_jitter, qtl=parkon_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))
#parkon_init_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_init_qtl)
#summary(fitqtl(parkon_jitter, qtl=parkon_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))

#parkon_expand_qtl<-addtoqtl(parkon_jitter,parkon_init_qtl, 2, 91)
#summary(fitqtl(parkon_jitter, qtl=parkon_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
#parkon_expand_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand_qtl)
#summary(fitqtl(parkon_jitter, qtl=parkon_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)) ## coalesces with first peak on chr 2


#parkon_expand_qtl<-addtoqtl(parkon_jitter,parkon_init_qtl, 5, 63)
#summary(fitqtl(parkon_jitter, qtl=parkon_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
#parkon_expand_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand_qtl)
#summary(fitqtl(parkon_jitter, qtl=parkon_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))

#parkon_expand_qtl_addint<-addint(parkon_jitter, pheno.col="pr",  qtl=parkon_expand_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)
#summary(fitqtl(parkon_jitter, qtl=parkon_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q2:Q3 + Q2:Q4 + Q3:Q5 + Q5:Q6))
#parkon_expand_qtl_add <- addqtl(parkon_jitter, qtl=parkon_expand_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6+ Q2:Q3 + Q2:Q4 + Q3:Q5 + Q5:Q6)

#parkon_expand2_qtl<-addtoqtl(parkon_jitter,parkon_expand_qtl, c(3,6), c(48.8,59.9))
#summary(fitqtl(parkon_jitter, qtl=parkon_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q2:Q3 + Q2:Q4 + Q3:Q5 + Q5:Q6))
#parkon_expand2_qtl<-refineqtl(parkon_jitter, pheno.col="pr", method="imp", qtl = parkon_expand2_qtl formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q2:Q3 + Q2:Q4 + Q3:Q5 + Q5:Q6)
#summary(fitqtl(parkon_jitter, qtl=parkon_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q2:Q3 + Q2:Q4 + Q3:Q5 + Q5:Q6))
 


## multiple QTL models / stepwiseqtl

#penalties=calc.penalties(twoqtl_parkon_hk_1000perm, alpha=0.05)
#stepwise_qtl_parkon_imp<-stepwiseqtl(parkon_jitter_sim, method="imp", model="normal", scan.pairs=FALSE, pheno.col="pr", penalties=penalties, keeptrace=TRUE)
#stepwise_qtl_parkon_hk<-stepwiseqtl(parkon_jitter_prob, method="hk", model="normal", scan.pairs=FALSE, pheno.col="pr", penalties=penalties, keeptrace=TRUE)


#stepwise_qtl_parkon_imp_dropQ367<-dropfromqtl(stepwise_qtl_parkon_imp, c(3,6,7))
#stepwise_qtl_parkon_imp_dropQ367<-refineqtl(parkon_jitter_sim,qtl=stepwise_qtl_parkon_imp_dropQ367,formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7, method="imp", pheno.col="pr")
#stepwise_qtl_parkon_imp_dropQ367_addint<-addint(parkon_jitter_sim, pheno.col="pr", qtl=stepwise_qtl_parkon_imp_dropQ367, method="imp")
#stepwise_qtl_parkon_imp_final<-refineqtl(parkon_jitter_sim, pheno.col="pr", method="imp", qtl = stepwise_qtl_parkon_imp_dropQ367, formula = y~Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q2:Q4+Q2:Q5+Q3:Q5)

#summary(fitqtl(parkon_jitter_sim, qtl=parkon_init_qtl, pheno.col="pr", method="imp"))
#summary(fitqtl(parkon_jitter_sim, qtl=stepwise_qtl_parkon_imp, formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q10 + Q2:Q8 + Q5:Q6 + Q1:Q4 + Q3:Q8 + Q7:Q9, pheno.col="pr", method="imp"))
#summary(fitqtl(parkon_jitter_sim, qtl=stepwise_qtl_parkon_imp_final, formula = y~Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q2:Q4+Q2:Q5+Q3:Q5,pheno.col="pr", method="imp"))


### prukoh ###
# prukoh has been coded the right way, L. kohalensis alleles are A, L. pruna alleles are B

prukoh<-read.cross(format="csv", file="PruKoh_Comprehensive_CrossFile.csv", estimate.map=FALSE)
prukoh_jitter<-jittermap(prukoh)
prukoh_jitter<-sim.geno(prukoh_jitter, step=1, n.draws=1000, error.prob=0.001)
prukoh_jitter<-calc.genoprob(prukoh_jitter, step=1, error.prob=0.001)

## scanone

qtl_prukoh_imp<-scanone(prukoh_jitter,method="imp", pheno.col="pr")
qtl_prukoh_hk<-scanone(prukoh_jitter,method="hk", pheno.col="pr")
qtl_prukoh_ehk<-scanone(prukoh_jitter,method="ehk", pheno.col="pr")
qtl_prukoh_em<-scanone(prukoh_jitter,method="em", pheno.col="pr")

qtl_prukoh_imp_1000perm<-scanone(prukoh_jitter,method="imp", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE) # this takes about 2 minutes
qtl_prukoh_hk_1000perm<-scanone(prukoh_jitter,method="hk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
qtl_prukoh_em_1000perm<-scanone(prukoh_jitter,method="em", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
#qtl_prukoh_ehk_1000perm<-scanone(prukoh_jitter,method="ehk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)


## cim 

#cim_prukoh_imp<-cim(prukoh_jitter, pheno.col="pr", method="imp", n.marcovar=3, window=10)
cim_prukoh_hk<-cim(prukoh_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20)
#cim_prukoh_ehk<-cim(prukoh_jitter, pheno.col="pr", method="ehk", n.marcovar=3, window=10)
#cim_prukoh_em<-cim(prukoh_jitter, pheno.col="pr", method="em", n.marcovar=3, window=10)

cim_prukoh_hk_1000perm<-cim(prukoh_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20, n.perm=1000)

cairo_pdf("cim_prukoh_hk_5_20.pdf", width=8, height=3.5)
plot(cim_prukoh_hk, ylim=c(0,35));abline(h=summary(cim_prukoh_hk_1000perm)[1,1])
dev.off()

qtl_prukoh<-c("1","2","3","4","5","6","X")

## scantwo
twoqtl_prukoh_imp<-scantwo(prukoh_jitter, pheno.col="pr", method="imp", chr=qtl_prukoh, clean.output=TRUE)
twoqtl_prukoh_hk<-scantwo(prukoh_jitter, pheno.col="pr", method="hk", chr=qtl_prukoh, clean.output=TRUE)

plot(twoqtl_prukoh_imp, lower="cond-int") # strong evidence for interaction chr1:chr2
plot(twoqtl_prukoh_hk, lower="cond-int")

plot(twoqtl_prukoh_hk, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_hk, chr=2, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_hk, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_hk, chr=4, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_hk, chr=5, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_hk, chr=6, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_hk, chr="X", lower="cond-int", upper="cond-add")

# permutate the scantwo function for LOD  penalties
twoqtl_prukoh_hk_1000perm<-scantwo(prukoh_jitter, pheno.col="pr", method="hk", chr=qtl_prukoh, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes
summary(twoqtl_prukoh_hk, threshold=as.matrix(as.data.frame(summary(twoqtl_prukoh_hk_10000perm, 0.05)[1:5])[1,]))
penalties=calc.penalties(twoqtl_prukoh_hk_1000perm, alpha=0.05)


    main    heavy    light 
3.389887 5.461719 3.136901

max(qtl_prukoh_hk)

prukoh_init_qtl<-makeqtl(prukoh_jitter,chr=2, pos=71)
summary(fitqtl(prukoh_jitter, qtl=prukoh_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1))
prukoh_init_qtl_add <- addqtl(prukoh_jitter, qtl=prukoh_init_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1)

prukoh_expand_qtl<-addtoqtl(prukoh_jitter,prukoh_init_qtl, 1, 46)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
prukoh_expand_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand_qtl)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
prukoh_expand_qtl_addint<-addint(prukoh_jitter, pheno.col="pr",  qtl=prukoh_expand_qtl, method="imp", formula=  y ~ Q1 + Q2)
prukoh_expand_qtl_add <- addqtl(prukoh_jitter, qtl=prukoh_expand_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2)


prukoh_expand2_qtl<-addtoqtl(prukoh_jitter,prukoh_expand_qtl, 6, 123)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
prukoh_expand2_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand2_qtl)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
prukoh_expand2_qtl_addint<-addint(prukoh_jitter, pheno.col="pr",  qtl=prukoh_expand2_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3)
prukoh_expand2_qtl_qtl_add <- addqtl(prukoh_jitter, qtl=prukoh_expand2_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3)

prukoh_expand3_qtl<-addtoqtl(prukoh_jitter,prukoh_expand2_qtl, 5, 87)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4))
prukoh_expand3_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand3_qtl)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4))
prukoh_expand3_qtl_addint<-addint(prukoh_jitter, pheno.col="pr",  qtl=prukoh_expand3_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4)
prukoh_expand3_qtl_add <- addqtl(prukoh_jitter, qtl=prukoh_expand3_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4)


prukoh_expand4_qtl<-addtoqtl(prukoh_jitter,prukoh_expand3_qtl, "X", 12)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))
prukoh_expand4_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand4_qtl)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))
prukoh_expand4_qtl_addint<-addint(prukoh_jitter, pheno.col="pr",  qtl=prukoh_expand4_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5)
prukoh_expand4_qtl_add <- addqtl(prukoh_jitter, qtl=prukoh_expand4_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5)

prukoh_expand5_qtl<-addtoqtl(prukoh_jitter,prukoh_expand4_qtl, 4, 85)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
prukoh_expand5_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand5_qtl)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
prukoh_expand5_qtl_addint<-addint(prukoh_jitter, pheno.col="pr",  qtl=prukoh_expand5_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)
prukoh_expand5_qtl_add <- addqtl(prukoh_jitter, qtl=prukoh_expand5_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)

prukoh_expand6_qtl<-addtoqtl(prukoh_jitter,prukoh_expand5_qtl, 3,79)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7))
prukoh_expand6_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand6_qtl)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7))
prukoh_expand6_qtl_addint<-addint(prukoh_jitter, pheno.col="pr",  qtl=prukoh_expand6_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7)
prukoh_expand6_qtl_add <- addqtl(prukoh_jitter, qtl=prukoh_expand6_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7)

prukoh_expand7_qtl<-addtoqtl(prukoh_jitter,prukoh_expand6_qtl, 2, 140)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8))
prukoh_expand7_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand7_qtl)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8))
prukoh_expand7_qtl_addint<-addint(prukoh_jitter, pheno.col="pr",  qtl=prukoh_expand7_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8)
prukoh_expand7_qtl_add <- addqtl(prukoh_jitter, qtl=prukoh_expand7_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8)

prukoh_expand8_qtl<-addtoqtl(prukoh_jitter,prukoh_expand7_qtl, 2, 28)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9))
prukoh_expand8_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand8_qtl)
summary(fitqtl(prukoh_jitter, qtl=prukoh_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9))

prukoh_QTLintervals<-data.frame(rbind(
cbind(paste("LODint1"),as.character(lodint(prukoh_expand7_qtl, qtl.index=1, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=1, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=1, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(prukoh_expand7_qtl, qtl.index=1, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=1, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=1, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(prukoh_expand7_qtl, qtl.index=2, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=2, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=2, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(prukoh_expand7_qtl, qtl.index=2, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=2, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=2, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(prukoh_expand7_qtl, qtl.index=3, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=3, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=3, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(prukoh_expand7_qtl, qtl.index=3, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=3, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=3, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(prukoh_expand7_qtl, qtl.index=4, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=4, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=4, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(prukoh_expand7_qtl, qtl.index=4, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=4, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=4, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(prukoh_expand7_qtl, qtl.index=5, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=5, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=5, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(prukoh_expand7_qtl, qtl.index=5, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=5, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=5, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(prukoh_expand7_qtl, qtl.index=6, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=6, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=6, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(prukoh_expand7_qtl, qtl.index=6, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=6, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=6, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(prukoh_expand7_qtl, qtl.index=7, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=7, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=7, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(prukoh_expand7_qtl, qtl.index=7, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=7, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=7, drop=2, expandtomarker=TRUE))))]),
cbind(paste("LODint1"),as.character(lodint(prukoh_expand7_qtl, qtl.index=8, drop=1, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=8, drop=1, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=8, drop=1, expandtomarker=TRUE))))]),
cbind(paste("LODint2"),as.character(lodint(prukoh_expand7_qtl, qtl.index=8, drop=2, expandtomarker=TRUE)[1,1]),rownames(lodint(prukoh_expand7_qtl, qtl.index=8, drop=2, expandtomarker=TRUE))[c(1,length(rownames(lodint(prukoh_expand7_qtl, qtl.index=8, drop=2, expandtomarker=TRUE))))]))
)


prukoh_newLG_QTLeffects<-data.frame(rbind(
data.frame(marker=find.marker(prukoh_newLG_jitter, 2,54),markerpos=find.markerpos(prukoh_newLG_jitter,as.character(find.marker(prukoh_newLG_jitter, 2,54)))$pos,matrix(effectplot(prukoh_newLG_jitter,"2@54.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(prukoh_newLG_jitter, 2,127.5),markerpos=find.markerpos(prukoh_newLG_jitter,as.character(find.marker(prukoh_newLG_jitter, 2,127.5)))$pos,matrix(effectplot(prukoh_newLG_jitter,"2@127.5", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(prukoh_newLG_jitter, 3,104),markerpos=find.markerpos(prukoh_newLG_jitter,as.character(find.marker(prukoh_newLG_jitter, 3,104)))$pos,matrix(effectplot(prukoh_newLG_jitter,"3@104.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(prukoh_newLG_jitter, 1,61.5),markerpos=find.markerpos(prukoh_newLG_jitter,as.character(find.marker(prukoh_newLG_jitter, 1,61.5)))$pos,matrix(effectplot(prukoh_newLG_jitter,"1@61.5", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(prukoh_newLG_jitter, 4,81),markerpos=find.markerpos(prukoh_newLG_jitter,as.character(find.marker(prukoh_newLG_jitter, 4,81)))$pos,matrix(effectplot(prukoh_newLG_jitter,"4@81.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(prukoh_newLG_jitter, 5,87),markerpos=find.markerpos(prukoh_newLG_jitter,as.character(find.marker(prukoh_newLG_jitter, 5,87)))$pos,matrix(effectplot(prukoh_newLG_jitter,"5@87.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(prukoh_newLG_jitter, 6,70.0),markerpos=find.markerpos(prukoh_newLG_jitter,as.character(find.marker(prukoh_newLG_jitter, 6,70.0)))$pos,matrix(effectplot(prukoh_newLG_jitter,"6@70.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)),
data.frame(marker=find.marker(prukoh_newLG_jitter, "X",14.0),markerpos=find.markerpos(prukoh_newLG_jitter,as.character(find.marker(prukoh_newLG_jitter, "X",14.0)))$pos,matrix(effectplot(prukoh_newLG_jitter,"X@14.0", pheno.col="pr", draw=FALSE)$Means,ncol=3)))
)


prukoh_QTLeffects$est_single<-(prukoh_QTLeffects$X1-prukoh_QTLeffects$X3)/2
prukoh_QTLeffects$est_single[8]<-(prukoh_QTLeffects$X1[8]-prukoh_QTLeffects$X2[8])/2
prukoh_QTLeffects$est<-summary(fitqtl(prukoh_jitter, qtl=prukoh_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8))$ests[c(2,15,13,4,11,8,6,10),1]*(-1)

prukoh_QTLeffects$PVE<-100*prukoh_QTLeffects$est/2
prukoh_QTLeffects$QTL<-c("2a","2b","3a","1","4","5b","6","X")

peakmarkerlist_prukoh=c(substr(as.matrix(find.flanking(prukoh_newLG_jitter,1,61.5)),0,7),
substr(as.matrix(find.flanking(prukoh_newLG_jitter,2,54.0)),0,7),
substr(as.matrix(find.flanking(prukoh_newLG_jitter,2,127.5)),0,7),
substr(as.matrix(find.flanking(prukoh_newLG_jitter,3,104.0)),0,7),
substr(as.matrix(find.flanking(prukoh_newLG_jitter,4,81.0)),0,7),
substr(as.matrix(find.flanking(prukoh_newLG_jitter,5,87.0)),0,7),
substr(as.matrix(find.flanking(prukoh_newLG_jitter,6,70.0)),0,7),
substr(as.matrix(find.flanking(prukoh_newLG_jitter,"X",14.0)),0,7))

write.table(peakmarkerlist_prukoh,"peakmarkerlist_prukoh.txt", quote=FALSE, row.names=FALSE, sep="\t")


prukoh_QTLbayesintervals<-data.frame(LG=character(), locus=character(), position=numeric())
CIsize_prukoh<-c()
for(i in 1:8) {
	temp_df<-data.frame(LG= bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1],locus=names(pull.map(prukoh_newLG_jitter)[[as.character(bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]),position=matrix(pull.map(prukoh_newLG_jitter)[[as.character(bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]))

	from=min(bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	to=max(bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	
	CIsize_prukoh=c(CIsize_prukoh,to-from)

	prukoh_QTLbayesintervals<-rbind(prukoh_QTLbayesintervals,temp_df[max(which(temp_df$position<=from)):min(which(temp_df$position>=to)),])
	}

prukoh_QTLbayesintervals$scaffold<-gsub("_.*","",prukoh_QTLbayesintervals$locus)

write.table(prukoh_QTLbayesintervals,"prukoh_QTLbayesintervals.txt", quote=FALSE, row.names=FALSE, sep="\t")

additive_effects_prukoh<-matrix(summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8))$ests[c(2,4,6,8,10,11,13,15),1])
additive_effects_prukoh<-(-1)*additive_effects_prukoh

summary(lm(CIsize_prukoh~additive_effects_prukoh))

all_additive_effects<-rbind(additive_effects_parkoh,additive_effects_parkon,additive_effects_prukoh,additive_effects_LPmerge)
all_CIsizes<-c(CIsize_parkoh,CIsize_parkon,CIsize_prukoh, CIsize_cereuk)
cross_list<-c(rep("parkoh",9),rep("parkon",9),rep("prukoh",8),rep("cereuk",7))

summary(lm(all_CIsizes~all_additive_effects + cross_list))

cairo_pdf("PruKoh_fitqtl.pdf", width=8, height=3.5)
plotLodProfile(prukoh_newLG_expand7_qtl,showallchr=TRUE, qtl.labels=FALSE); abline(h=3.38)
dev.off()

fitdistr(prukoh_QTLeffects$est, "exponential") # fit exponential distribution to data
cairo_pdf("prukoh_effectsize_distr.pdf", width=3, height=3)
plot(density(rexp(10000, rate=fitdistr(prukoh_QTLeffects$est, "exponential")$estimate)),xlim=c(0,1))
hist(prukoh_QTLeffects$est, prob=FALSE, add=TRUE, breaks=4)
dev.off()


fitdistr(cereuk_QTLeffects$est, "gamma") # fit gamma distribution to data
cairo_pdf("cereuk_effectsize_distr.pdf", width=3, height=3)
plot(density(rgamma(10000,shape=3.12, rate=37.13)),xlim=c(0,1))
hist(cereuk_QTLeffects$est, prob=TRUE, breaks=4, add=TRUE)
dev.off()

true_loci_prukoh<-otto_jones_ci(D=0.975,M=0.091081925,nd=8,alpha=0.05,amin=0.0479672,res=4,res2=2,max.loci=100)

######### probabilities of observed QTL overlap

qtl_overlap<-function(n_qtl=15,n_1=8,n_2=9,n_shared=6,res=1e5) {
	overlapping<-c()
	for(i in 1:res) {
		qtl1<-sample(c(1:n_qtl),n_1)
		qtl2<-sample(c(1:n_qtl),n_2)
		overlapping<-c(overlapping,length(which(qtl1 %in% qtl2)))
				}
	probability=(which(quantile(overlapping, seq(0,1,1/res))>=n_shared)[1]-1)/res
	return(1-probability)		
	}

qtl_overlap3<-function(n_qtl=12,n_1=8,n_2=9, n_3=7, n_shared=5, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) {
	overlapping<-c()
	unique_1<-c()
	unique_2<-c()
	unique_3<-c()
	for(i in 1:res) {
		qtl1<-sample(c(1:n_qtl),n_1)
		qtl2<-sample(c(1:n_qtl),n_2)
		qtl3<-sample(c(1:n_qtl),n_3)
		share12<-qtl1[which(qtl1 %in% qtl2)]
		overlapping<-c(overlapping,length(which(share12 %in% qtl3)))
		unique_1<-c(unique_1,length(qtl1)-length(which(qtl1 %in% qtl2 | qtl1 %in% qtl3)))
		unique_2<-c(unique_2,length(qtl2)-length(which(qtl2 %in% qtl1 | qtl2 %in% qtl3)))
		unique_3<-c(unique_3,length(qtl3)-length(which(qtl3 %in% qtl1 | qtl3 %in% qtl2)))
		}
	probability_sharing=1-(which(quantile(overlapping, seq(0,1,1/res))>=n_shared)[1]-1)/res
	probability_unique1=(max(which(quantile(unique_1, seq(0,1,1/res))<=n_unique1))-1)/res
	probability_unique2=(max(which(quantile(unique_2, seq(0,1,1/res))<=n_unique2))-1)/res
	probability_unique3=(max(which(quantile(unique_3, seq(0,1,1/res))<=n_unique3))-1)/res
	out=list(probability_sharing,probability_unique1,probability_unique2,probability_unique3)
	return(out)		
	}

# average ParKon n_qtl estimate, 12 QTL:
qtl_overlap(n_qtl=12,n_1=9,n_2=8,n_shared=6,res=1e5)
#P = 0.74
qtl_overlap(n_qtl=12,n_1=9,n_2=7,n_shared=6,res=1e5)
#P = 0.36
qtl_overlap(n_qtl=12,n_1=7,n_2=8,n_shared=6,res=1e5)
#P = 0.15

# average ParKon n_qtl estimate, 16 QTL:
qtl_overlap(n_qtl=16,n_1=9,n_2=8,n_shared=6,res=1e5)
#P = 0.
qtl_overlap(n_qtl=16,n_1=9,n_2=7,n_shared=6,res=1e5)
#P = 0.
qtl_overlap(n_qtl=16,n_1=7,n_2=8,n_shared=6,res=1e5)
#P = 0.

# average PruKoh n_qtl estimate, 21 QTL:
qtl_overlap(n_qtl=21,n_1=9,n_2=8,n_shared=6,res=1e5)
#P = 0.028
qtl_overlap(n_qtl=21,n_1=9,n_2=7,n_shared=6,res=1e5)
#P = 0.009
qtl_overlap(n_qtl=21,n_1=7,n_2=8,n_shared=6,res=1e5)
#P = 0.004


# minimum triple and unique:
qtl_overlap3(n_qtl=12,n_1=8,n_2=9, n_3=7, n_shared=5, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 


qtl_overlap3(n_qtl=12,n_1=8,n_2=9, n_3=7, n_shared=4, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 

# mean triple and unique:
qtl_overlap3(n_qtl=16,n_1=8,n_2=9, n_3=7, n_shared=5, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 


qtl_overlap3(n_qtl=16,n_1=8,n_2=9, n_3=7, n_shared=4, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 

# maximum triple and unique:
qtl_overlap3(n_qtl=21,n_1=8,n_2=9, n_3=7, n_shared=5, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 


qtl_overlap3(n_qtl=21,n_1=8,n_2=9, n_3=7, n_shared=4, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 


# average CerEuk estimate, 15
#qtl_overlap(n_qtl=15,n_1=9,n_2=7,n_shared=6,res=1e5)
#P = 0.084	

#qtl_overlap(n_qtl=15,n_1=7,n_2=8,n_shared=6,res=1e5)
#P = 0.032


## overlap in the sense of randomly distributed windows:

average_windowsize=(11+29+41+19+40+26+70+8+14+8+46+14+1+15+2+4)/16 # based on ParKon and PruKoh QTL, 95 Bayesian CI

LG_length_parkon<-c();for(i in 1:8) { LG_length_parkon<-c(LG_length_parkon,max(pull.map(parkon)[[i]]))}; sum(LG_length_parkon) # total map length ParKon
LG_length_prukoh<-c();for(i in 1:8) { LG_length_prukoh<-c(LG_length_prukoh,max(pull.map(prukoh)[[i]]))}; sum(LG_length_prukoh) # total map length PruKoh

average_LG_length=((sum(LG_length_parkon)/8)+(sum(LG_length_prukoh)/8))/2

qtl_overlap_randomwindows<-function(chr_size=90,window_size=22,qtl_number=9,shared_qtl=6,nboot=100000) {
	overlapping<-c()
	for(n in 1:nboot) {
		sample_1_start<-sample(1:chr_size,1)
		sample_2_start<-sample(1:chr_size,1)
		sample1<-seq(1:chr_size)[sample_1_start:(sample_1_start+window_size)]
		sample2<-seq(1:chr_size)[sample_2_start:(sample_2_start+window_size)]
		if(length(which(sample1 %in% sample2)) > 1) { overlapping<-c(overlapping,n) }
		}
	single_prob=length(overlapping)/nboot
	iterative_prob=(single_prob^shared_qtl)*((1-single_prob)^(qtl_number-shared_qtl)) * ncol(combn(qtl_number, shared_qtl))
	output<-list()
	output$iterative_prob<-iterative_prob
	output$overlap_prob<-single_prob
	return(output)
	}

qtl_overlap_randomwindows(chr_size=average_LG_length,window_size=average_windowsize,qtl_number=9,shared_qtl=6,nboot=100000)
0.02463777
qtl_overlap_randomwindows(chr_size=average_LG_length,window_size=average_windowsize,qtl_number=8,shared_qtl=6,nboot=100000)
0.01198767
qtl_overlap_randomwindows(chr_size=average_LG_length,window_size=average_windowsize,qtl_number=7,shared_qtl=6,nboot=100000)
0.004445523


qtl_overlap_randomwindows(chr_size=150,window_size=10,qtl_number=8,shared_qtl=6,nboot=100)

# prukoh_init_qtl<-makeqtl(prukoh_jitter, chr=c(1,2,5,6), pos=c(48,71,87.9,100))
# summary(fitqtl(prukoh_jitter, qtl=prukoh_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4))
# prukoh_init_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_init_qtl)
# summary(fitqtl(prukoh_jitter, qtl=prukoh_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4))
# summary(fitqtl(prukoh_jitter, qtl=prukoh_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q1:Q2))

#prukoh_expand_qtl<-addtoqtl(prukoh_jitter,prukoh_init_qtl, 2, 113)
#summary(fitqtl(prukoh_jitter, qtl=prukoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
#prukoh_expand_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand_qtl)
#summary(fitqtl(prukoh_jitter, qtl=prukoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))


# prukoh_add_qtl <- addqtl(prukoh_jitter, qtl=prukoh_expand_qtl, pheno.col="pr", method="imp")
# prukoh_expand_qtl<-addtoqtl(prukoh_jitter,prukoh_init_qtl, c(3,4,"X"), c(100, 83, 12))
# prukoh_expand_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand_qtl)
# summary(fitqtl(prukoh_jitter, qtl=prukoh_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q2))

# prukoh_expand_qtl_addint<-addint(prukoh_jitter, pheno.col="pr",  qtl=prukoh_expand_qtl, method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7)
# prukoh_add_qtl <- addqtl(prukoh_jitter, qtl=prukoh_expand_qtl, pheno.col="pr", method="imp",formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7)

# prukoh_expand2_qtl<-addtoqtl(prukoh_jitter,prukoh_expand_qtl, 2, 140)
# prukoh_expand2_qtl<-refineqtl(prukoh_jitter, pheno.col="pr", method="imp", qtl = prukoh_expand2_qtl)
# summary(fitqtl(prukoh_jitter, qtl=prukoh_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8)

# prukoh_expand2_qtl_addint<-addint(prukoh_jitter, pheno.col="pr",  qtl=prukoh_expand2_qtl, method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 +Q8)

# summary(fitqtl(prukoh_jitter, qtl=prukoh_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8))

# multiple QTL models / stepwiseqtl
# stepwise_qtl_prukoh_imp<-stepwiseqtl(prukoh_jitter_sim, method="imp", model="normal", scan.pairs=FALSE, pheno.col="pr", penalties=penalties, keeptrace=TRUE)
# stepwise_qtl_prukoh_hk<-stepwiseqtl(prukoh_jitter_prob, method="hk", model="normal", scan.pairs=FALSE, pheno.col="pr", penalties=penalties, keeptrace=TRUE)

# stepwise_qtl_prukoh_imp_dropQ68<-dropfromqtl(stepwise_qtl_prukoh_imp, c(6,8))
# stepwise_qtl_prukoh_imp_dropQ68<-refineqtl(prukoh_jitter_sim,qtl=stepwise_qtl_prukoh_imp_dropQ68,formula=y~Q1+Q2+Q3+Q4+Q5+Q7+Q8, method="imp", pheno.col="pr")
# stepwise_qtl_prukoh_imp_dropQ68_addint<-addint(prukoh_jitter_sim, pheno.col="pr", qtl=stepwise_qtl_prukoh_imp_dropQ68, method="imp")
# stepwise_qtl_prukoh_imp_final<-refineqtl(prukoh_jitter_sim, pheno.col="pr", method="imp", qtl = stepwise_qtl_prukoh_imp_dropQ68, formula = y~Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q1:Q2)

# summary(fitqtl(prukoh_jitter_sim, qtl=prukoh_init_qtl, pheno.col="pr", method="imp"))
# summary(fitqtl(prukoh_jitter_sim, qtl=stepwise_qtl_prukoh_imp, formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q10 + Q4:Q5 + Q2:Q9 + Q1:Q6, pheno.col="pr", method="imp"))
# summary(fitqtl(prukoh_jitter_sim, qtl=stepwise_qtl_prukoh_imp_dropQ68, formula = y~Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8,pheno.col="pr", method="imp"))
# summary(fitqtl(prukoh_jitter_sim, qtl=stepwise_qtl_prukoh_imp_final, formula = y~Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q1:Q2,pheno.col="pr", method="imp"))


############################### match lodprofile peaks

load("CerEuk_LPmerge.RData")
load("CerEuk_LPmerge_expand6_qtl.RData")
load("cereuk_QTLeffects.RData")

load("CerEuk_HC1.RData")
load("CerEukHC1_expand2_qtl.RData")

# combinedmap<-data.frame(group=factor(),position=numeric(),locus=character())
# for(i in 1:length(pull.map(parkon))) { 
	# combinedmap=rbind(combinedmap,data.frame(group=paste0("parkonLG",i),position=matrix(pull.map(parkon)[[i]]),locus=substr(names(pull.map(parkon)[[i]]),0,7)))
	# }
# for(i in 1:length(pull.map(prukoh))) { 
	# combinedmap=rbind(combinedmap,data.frame(group=paste0("prukohLG",i),position=matrix(pull.map(prukoh)[[i]]),locus=substr(names(pull.map(prukoh)[[i]]),0,7)))
	# }
# #for(i in 1:length(pull.map(parkoh))) { 
# #	combinedmap=rbind(combinedmap,data.frame(group=paste0("parkohLG",i),position=matrix(pull.map(parkoh)[[i]]),locus=substr(names(pull.map(parkoh)[[i]]),0,7)))
# #	}
# for(i in 1:length(pull.map(CerEuk_LPmerge))) { 
	# combinedmap=rbind(combinedmap,data.frame(group=paste0("cereukLG",i),position=matrix(pull.map(CerEuk_LPmerge)[[i]]),locus=substr(names(pull.map(CerEuk_LPmerge)[[i]]),0,7)))
	# }

# load r function make_mapchart from Laupala Stuff/Scripts/make_mapchart_from_rQTL.r
# make MapChart files for Bayes CI comparisons	
mapthese=c("parkonLG1","prukohLG1","parkonLG2","prukohLG2","parkonLG3","prukohLG3","parkonLG4","prukohLG4","parkonLG5","prukohLG5","parkonLG6","prukohLG6","parkonLG7","prukohLG7","parkonLG8","prukohLG8")
make_mapchart(cross_objects=list(parkon,prukoh), name=c("parkon","prukoh"), qtl_results=list(parkon_expand8_qtl,prukoh_expand7_qtl), mapthese=mapthese,output="parkon_prukoh_combinedQTL.mct")

mapthese=c("parkonLG1","cereukLG1","parkonLG2","cereukLG2","parkonLG3","cereukLG3","parkonLG4","cereukLG4","parkonLG5","cereukLG5","parkonLG6","cereukLG6","parkonLG7","cereukLG7","parkonLG8","cereukLG8")
make_mapchart(cross_objects=list(parkon,CerEuk_LPmerge), name=c("parkon","cereuk"), qtl_results=list(parkon_expand8_qtl,CerEuk_LPmerge_expand6_qtl), mapthese=mapthese,output="parkon_cereuk_combinedQTL.mct")

mapthese=c("cereukLG1","prukohLG1","cereukLG2","prukohLG2","cereukLG3","prukohLG3","cereukLG4","prukohLG4","cereukLG5","prukohLG5","cereukLG6","prukohLG6","cereukLG7","prukohLG7","cereukLG8","prukohLG8")
make_mapchart(cross_objects=list(CerEuk_LPmerge,prukoh), name=c("cereuk","prukoh"), qtl_results=list(CerEuk_LPmerge_expand6_qtl,prukoh_expand7_qtl), mapthese=mapthese,output="cereuk_prukoh_combinedQTL.mct")

# library(LinkageMapView)

# outfile=file.path(getwd(),"parkon_prukoh_QTLoverlap_BayInt.pdf")
# outfile=file.path(getwd(),"parkoh_cereuk_QTLoverlap_BayInt.pdf")

# peakmarkerlist_parkon=c(substr(as.matrix(find.flanking(parkon_jitter,1,79.1)),0,7),
# substr(as.matrix(find.flanking(parkon_jitter,2,57.6)),0,7),
# substr(as.matrix(find.flanking(parkon_jitter,3,60.7)),0,7),
# substr(as.matrix(find.flanking(parkon_jitter,4,48.0)),0,7),
# substr(as.matrix(find.flanking(parkon_jitter,5,32.5)),0,7),
# substr(as.matrix(find.flanking(parkon_jitter,5,59.0)),0,7),
# substr(as.matrix(find.flanking(parkon_jitter,"X",60.0)),0,7))

# peakmarkerlist_prukoh=c(substr(as.matrix(find.flanking(prukoh_jitter,1,61.0)),0,7),
# substr(as.matrix(find.flanking(prukoh_jitter,2,71.0)),0,7),
# substr(as.matrix(find.flanking(prukoh_jitter,2,140)),0,7),
# substr(as.matrix(find.flanking(prukoh_jitter,3,79.0)),0,7),
# substr(as.matrix(find.flanking(prukoh_jitter,4,83.0)),0,7),
# substr(as.matrix(find.flanking(prukoh_jitter,5,88.0)),0,7),
# substr(as.matrix(find.flanking(prukoh_jitter,6,69.0)),0,7),
# substr(as.matrix(find.flanking(prukoh_jitter,"X",13.0)),0,7))


# peakmarkerlist_parkoh=c(substr(as.matrix(find.flanking(parkoh_jitter,2,70.0)),0,7),
# substr(as.matrix(find.flanking(parkoh_jitter,2,198.1)),0,7),
# substr(as.matrix(find.flanking(parkoh_jitter,3,74.0)),0,7),
# substr(as.matrix(find.flanking(parkoh_jitter,3,165.6)),0,7),
# substr(as.matrix(find.flanking(parkoh_jitter,1,86.3)),0,7),
# substr(as.matrix(find.flanking(parkoh_jitter,4,39.0)),0,7),
# substr(as.matrix(find.flanking(parkoh_jitter,5,31.0)),0,7),
# substr(as.matrix(find.flanking(parkoh_jitter,6,57.0)),0,7),
# substr(as.matrix(find.flanking(parkoh_jitter,"X",41.0)),0,7))

# peakmarkerlist_cereuk=c(substr(as.matrix(find.flanking(CerEuk_LPmerge_jitter,1,79.9)),0,7),
# substr(as.matrix(find.flanking(CerEuk_LPmerge_jitter,2,59)),0,7),
# substr(as.matrix(find.flanking(CerEuk_LPmerge_jitter,3,71)),0,7),
# substr(as.matrix(find.flanking(CerEuk_LPmerge_jitter,4,59.8)),0,7),
# substr(as.matrix(find.flanking(CerEuk_LPmerge_jitter,5,36.9)),0,7),
# substr(as.matrix(find.flanking(CerEuk_LPmerge_jitter,6,5)),0,7),
# substr(as.matrix(find.flanking(CerEuk_LPmerge_jitter,"X",38)),0,7))

# peakmarkerlist_cereuk[5]<-"S001891" # LG 2 has no markers in common with 

# sharedmarkers<-as.vector(combinedmap[substr(combinedmap[,"group"],0,6)=="parkon",][which(combinedmap[substr(combinedmap[,"group"],0,6)=="parkon","locus"] %in% combinedmap[substr(combinedmap[,"group"],0,6)=="prukoh","locus"]),"locus"])
# sharedmarkers<-as.vector(combinedmap[substr(combinedmap[,"group"],0,6)=="parkoh",][which(combinedmap[substr(combinedmap[,"group"],0,6)=="parkoh","locus"] %in% combinedmap[substr(combinedmap[,"group"],0,6)=="cereuk","locus"]),"locus"])


#lmv.linkage.plot(combinedmap, outfile, autoconnadj=TRUE, mapthese=c("parkonLG1","prukohLG1","parkonLG2","prukohLG2"),dupnbr=TRUE, showonly=peakmarkerlist)

# # qtldf<-data.frame(chr=character(), qtl=character(), so=numeric(),si=numeric(),ei=numeric(),eo=numeric(), col=character())

# # for (i in 1:9) {
	# # qtlone <- bayesint(reorderqtl(parkon_expand8_qtl), qtl.index = i)
	# # qtldf <- rbind(qtldf, data.frame( chr = paste0("parkonLG",as.character(qtlone[1,1])), qtl = paste0(qtlone[1,1], "@" , round(qtlone[2,2], digits = 1)),
        # # so = qtlone[1, 2],
        # # si = qtlone[2, 2],
        # # ei = qtlone[2, 2],
        # # eo = qtlone[3, 2],
        # # col="black"))
	# # }
# # for (i in 1:8) {
	# # qtlone <- bayesint(reorderqtl(prukoh_expand7_qtl), qtl.index = i)
	# # qtldf <- rbind(qtldf, data.frame( chr = paste0("prukohLG",as.character(qtlone[1,1])), qtl = paste0(qtlone[1,1], "@" , round(qtlone[2,2], digits = 1)),
        # # so = qtlone[1, 2],
        # # si = qtlone[2, 2],
        # # ei = qtlone[2, 2],
        # # eo = qtlone[3, 2],
        # # col="black"))
	# # }
	
# # qtldf<-data.frame(chr=character(), qtl=character(), so=numeric(),si=numeric(),ei=numeric(),eo=numeric(), col=character())

# # for (i in 1:9) {
	# # qtlone <- bayesint(reorderqtl(parkon_expand8_qtl), qtl.index = i)
	# # qtldf <- rbind(qtldf, data.frame( chr = paste0("parkohLG",as.character(qtlone[1,1])), qtl = paste0(qtlone[1,1], "@" , round(qtlone[2,2], digits = 1)),
        # # so = qtlone[1, 2],
        # # si = qtlone[2, 2],
        # # ei = qtlone[2, 2],
        # # eo = qtlone[3, 2],
        # # col="black"))
	# # }
# # for (i in 1:7) {
	# # qtlone <- bayesint(reorderqtl(CerEuk_LPmerge_expand6_qtl), qtl.index = i)
	# # qtldf <- rbind(qtldf, data.frame( chr = paste0("cereukLG",as.character(qtlone[1,1])), qtl = paste0(qtlone[1,1], "@" , round(qtlone[2,2], digits = 1)),
        # # so = qtlone[1, 2],
        # # si = qtlone[2, 2],
        # # ei = qtlone[2, 2],
        # # eo = qtlone[3, 2],
        # # col="black"))
	# # }

# # qtldf<-data.frame(chr=character(), qtl=character(), so=numeric(),si=numeric(),ei=numeric(),eo=numeric(), col=character())

# # for (i in 1:9) {
	# # qtlone <- bayesint(reorderqtl(parkon_expand8_qtl), qtl.index = i)
	# # qtldf <- rbind(qtldf, data.frame( chr = paste0("parkohLG",as.character(qtlone[1,1])), qtl = paste0(qtlone[1,1], "@" , round(qtlone[2,2], digits = 1)),
        # # so = qtlone[1, 2],
        # # si = qtlone[2, 2],
        # # ei = qtlone[2, 2],
        # # eo = qtlone[3, 2],
        # # col="black"))
	# # }
# # for (i in 1:7) {
	# # qtlone <- bayesint(reorderqtl(CerEuk_LPmerge_expand6_qtl), qtl.index = i)
	# # qtldf <- rbind(qtldf, data.frame( chr = paste0("cereukLG",as.character(qtlone[1,1])), qtl = paste0(qtlone[1,1], "@" , round(qtlone[2,2], digits = 1)),
        # # so = qtlone[1, 2],
        # # si = qtlone[2, 2],
        # # ei = qtlone[2, 2],
        # # eo = qtlone[3, 2],
        # # col="black"))
	# # }

# # qtldf$chr=gsub("X","8",qtldf$chr)

# lmv.linkage.plot(combinedmap, outfile, autoconnadj=TRUE, mapthese=c("parkonLG1","prukohLG1","parkonLG2","prukohLG2"),dupnbr=TRUE, showonly=peakmarkerlist,qtldf=qtldf)
# # lmv.linkage.plot(combinedmap, outfile, autoconnadj=TRUE, mapthese=c("parkonLG1","prukohLG1","parkonLG2","prukohLG2","parkonLG3","prukohLG3","parkonLG4","prukohLG4",
# # "parkonLG5","prukohLG5","parkonLG6","prukohLG6","parkonLG7","prukohLG7","parkonLG8","prukohLG8"), showonly=sharedmarkers,dupnbr=TRUE,qtldf=qtldf)

# # lmv.linkage.plot(combinedmap, outfile, autoconnadj=TRUE, mapthese=c("parkohLG1","cereukLG1","parkohLG2","cereukLG2","parkohLG3","cereukLG3","parkohLG4","cereukLG4",
# # "parkohLG5","cereukLG5","parkohLG6","cereukLG6","parkohLG7","cereukLG7","parkohLG8","cereukLG8"),dupnbr=TRUE,qtldf=qtldf)


# # or make output for mapchart

# map1=c("parkonLG1","parkonLG2","parkonLG3","parkonLG4","parkonLG5","parkonLG6","parkonLG7","parkonLG8")
# map2=c("prukohLG1","prukohLG2","prukohLG3","prukohLG4","prukohLG5","prukohLG6","prukohLG7","prukohLG8")
# # map=c("parkonLG1","prukohLG1","parkonLG2","prukohLG2","parkonLG3","prukohLG3","parkonLG4","prukohLG4","parkonLG5","prukohLG5","parkonLG6","prukohLG6","parkonLG7","prukohLG7","parkonLG8","prukohLG8")
# output="parkon_prukoh_combinedQTL.mct"
# qtls=qtldf

# # make_mapchart<-function(map=map,qtls=qtldf,output="parkon_prukoh_combinedQTL.mct") {
# # write(paste0("; this automatically generated MapChart file contains the following linkage groups: ",map),file=output)
# # write(paste(" "),file=output,append=TRUE)
# # for(LG in map) {
	
	# # write(paste0("group ",LG),file=output,append=TRUE)
	# # write.table(cbind(as.character(combinedmap[which(combinedmap$group == LG),"locus"]),combinedmap[which(combinedmap$group == LG),"position"]),file=output,append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")
	# # if(nrow(qtls[which(qtls$chr==LG),c("qtl","so","si","ei","eo")]) > 0) {
		# # write(paste("qtls"),file=output,append=TRUE)
		# # write.table(qtls[which(qtls$chr==LG),c("qtl","so","si","ei","eo")],file=output,append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE,sep=" ")
		# # }
	# # write(paste(" "),file=output,append=TRUE)
	
	# # }

# # }

############################# match qtl effect sizes
library(ggplot2)

QTLeffects_sharing<-read.delim("6sppQTL_effectsizes_sharing.txt")

#Fig 3A:
# cross factor order: cereuk konpar parkoh prukoh
pps<-c(1.66,1.46,3.01,1.64)
detectedQTLnumber<-c(7,9,8,8)
totalQTLnumber<-c(17,13,14,20)
avg_effect<-aggregate(na.omit(QTLeffects_sharing)[,"effect"], by=list(na.omit(QTLeffects_sharing)[,"cross"]), FUN=mean)$x
lm_phenodist_effect<-lm(avg_effect~pps)

Residuals:
        1         2         3         4 
-0.001028 -0.005206 -0.000669  0.006903 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)  
(Intercept)  0.01696    0.01014   1.673   0.2363  
pps          0.04099    0.00497   8.247   0.0144 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.006175 on 2 degrees of freedom
Multiple R-squared:  0.9714,    Adjusted R-squared:  0.9572 
F-statistic: 68.01 on 1 and 2 DF,  p-value: 0.01439

Response: avg_effect
          Df     Sum Sq    Mean Sq F value  Pr(>F)  
pps        1 0.00259336 0.00259336  68.014 0.01439 *
Residuals  2 0.00007626 0.00003813 

phenodist_effect_plot<-ggplot(data.frame(avg_effect=avg_effect,pps=pps),aes(x=avg_effect,y=pps)) +
geom_point() +
geom_smooth(method="lm",fullrange = TRUE) +
scale_x_continuous(limits=c(0,0.15)) + 
scale_y_continuous(limits=c(0,3.5)) + 
theme(legend.position="none")

lm_phenodist_detectedQTLnumber<-lm(detectedQTLnumber~pps)

Residuals:
       1        2        3        4 
-1.03660  0.93749  0.13831 -0.03919 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)  
(Intercept)   8.2517     1.6308   5.060   0.0369 *
pps          -0.1296     0.7996  -0.162   0.8862  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.9935 on 2 degrees of freedom
Multiple R-squared:  0.01296,   Adjusted R-squared:  -0.4806 
F-statistic: 0.02625 on 1 and 2 DF,  p-value: 0.8862

lm_phenodist_totalQTLnumber<-lm(totalQTLnumber~pps)

Residuals:
      1       2       3       4 
 0.6010 -3.6814 -0.4925  3.5728 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)  
(Intercept)   18.743      6.022   3.112   0.0896 .
pps           -1.412      2.953  -0.478   0.6797  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 3.669 on 2 degrees of freedom
Multiple R-squared:  0.1026,    Adjusted R-squared:  -0.3461 
F-statistic: 0.2287 on 1 and 2 DF,  p-value: 0.6797

# Fig 3B:

QTLeffects_cereuk_prukoh<-na.omit(merge(QTLeffects_sharing[which(QTLeffects_sharing$cross=="cereuk"),],QTLeffects_sharing[which(QTLeffects_sharing$cross=="prukoh"),],"QTL",sort=FALSE))
colnames(QTLeffects_cereuk_prukoh)[c(4,7)]<-c("effect_cereuk","effect_prukoh")
QTLeffects_cereuk_konpar<-na.omit(merge(QTLeffects_sharing[which(QTLeffects_sharing$cross=="cereuk"),],QTLeffects_sharing[which(QTLeffects_sharing$cross=="konpar"),],"QTL",sort=FALSE))
colnames(QTLeffects_cereuk_konpar)[c(4,7)]<-c("effect_cereuk","effect_konpar")
QTLeffects_konpar_prukoh<-na.omit(merge(QTLeffects_sharing[which(QTLeffects_sharing$cross=="konpar"),],QTLeffects_sharing[which(QTLeffects_sharing$cross=="prukoh"),],"QTL",sort=FALSE))
colnames(QTLeffects_konpar_prukoh)[c(4,7)]<-c("effect_konpar","effect_prukoh")

c(0,0.225)
QTL_effects_plot<-ggplot() +
 geom_point(data=QTLeffects_cereuk_konpar,aes(x=effect_cereuk,y=effect_konpar), color="blue") +
 geom_smooth(data=QTLeffects_cereuk_konpar,aes(x=effect_cereuk,y=effect_konpar), method="lm",color="blue", fullrange=TRUE) +
 geom_point(data=QTLeffects_konpar_prukoh,aes(x=effect_konpar,y=effect_prukoh), color="red") +
 geom_smooth(data=QTLeffects_konpar_prukoh,aes(x=effect_konpar,y=effect_prukoh), method="lm",color="red", fullrange=TRUE) + 
 geom_point(data=QTLeffects_cereuk_prukoh,aes(x=effect_cereuk,y=effect_prukoh), color="black") +
 geom_smooth(data=QTLeffects_cereuk_prukoh,aes(x=effect_cereuk,y=effect_prukoh), method="lm",color="black", fullrange=TRUE) +
 xlab("effect size (pps) cereuk (black/blue) / konpar (red)") +
 ylab("effect size (pps) prukoh (black/red) / konpar (blue)") +
 theme(legend.position="none")

# Fig 3 C:
QTLeffects_sharing_plot<-ggplot(QTLeffects_sharing) + 
geom_point(aes(x=factor(sharing,levels=c("unique","double","triple")), y=effect, color=as.factor(cross))) +
#geom_smooth(data=data.frame(na.omit(QTLeffects_sharing)),aes(x=sharing, y=effect), method="lm") +
theme(legend.position="none")

abline(lm(effect~sharing, data=QTLeffects_sharing[which(QTLeffects_sharing$cross != "parkoh"),])
lm_size_sharing<-lm(effect~cross+sharing, data=QTLeffects_sharing[which(QTLeffects_sharing$cross != "parkoh"),])


Residuals:
      Min        1Q    Median        3Q       Max 
-0.068642 -0.034331 -0.009349  0.015487  0.149176 

Coefficients:
               Estimate Std. Error t value Pr(>|t|)  
(Intercept)    0.074914   0.028328   2.645    0.016 *
crosskonpar   -0.006880   0.029357  -0.234    0.817  
crossprukoh    0.008244   0.030042   0.274    0.787  
sharingtriple  0.023686   0.026510   0.893    0.383  
sharingunique -0.031348   0.035623  -0.880    0.390  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05788 on 19 degrees of freedom
  (12 observations deleted due to missingness)
Multiple R-squared:  0.1488,    Adjusted R-squared:  -0.03036 
F-statistic: 0.8306 on 4 and 19 DF,  p-value: 0.5222

Response: effect
          Df   Sum Sq   Mean Sq F value Pr(>F)
cross      2 0.001659 0.0008296  0.2476 0.7832
sharing    2 0.009472 0.0047362  1.4135 0.2677
Residuals 19 0.063662 0.0033506 

library(gridExtra)
pdf("6sppQTL_Fig3_QTLeffects_sharing.pdf", width=4, height=9.5)
grid.arrange(phenodist_effect_plot,QTL_effects_plot,QTLeffects_sharing_plot, ncol=1)
dev.off()

library(lme4)
lme_size_sharing<-lmer(effect~sharing + (1|cross), data=QTLeffects_sharing[which(QTLeffects_sharing$cross != "parkoh"),])
lme_null<-lmer(effect~ (1|cross), data=QTLeffects_sharing[which(QTLeffects_sharing$cross != "parkoh"),])
anova(lme_size_sharing,lme_null)

Models:
lme_null: effect ~ (1 | cross)
lme_size_sharing: effect ~ sharing + (1 | cross)
                 Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
lme_null          3 -64.397 -60.863 35.198  -70.397                         
lme_size_sharing  5 -63.906 -58.016 36.953  -73.906 3.5092      2      0.173


# parkoh_parkon_QTLeffects<-rbind(parkoh_QTLeffects[c(1,3,5,6,7,9),c("est","PVE","QTL")],parkon_QTLeffects[c(1,2,3,4,5,7),c("est","PVE","QTL")])

# anova(lm(matrix(parkoh_parkon_QTLeffects$est,ncol=2)[,1]~matrix(parkoh_parkon_QTLeffects$est,ncol=2)[,2])) 
# Response: matrix(parkoh_parkon_QTLeffects$est, ncol = 2)[, 1]
                                                    # Df   Sum Sq  Mean Sq F value  Pr(>F)  
# matrix(parkoh_parkon_QTLeffects$est, ncol = 2)[, 2]  1 0.057065 0.057065   16.66 0.01508 *

# Coefficients
# (Intercept)  matrix(parkoh_parkon_QTLeffects$est, ncol = 2)[, 2]  
# -0.005293                                             1.531640

# anova(lm(matrix(parkoh_parkon_QTLeffects$PVE,ncol=2)[,1]~matrix(parkoh_parkon_QTLeffects$PVE,ncol=2)[,2])) 


# parkoh_prukoh_QTLeffects<-rbind(parkoh_QTLeffects[c(1,2,3,5,6,8,9),c("est","PVE","QTL")],prukoh_QTLeffects[c(1,2,3,4,5,7,8),c("est","PVE","QTL")])
# anova(lm(matrix(parkoh_prukoh_QTLeffects$est,ncol=2)[,1]~matrix(parkoh_prukoh_QTLeffects$est,ncol=2)[,2])) 
# Response: matrix(parkoh_prukoh_QTLeffects$est, ncol = 2)[, 1]
                                                    # Df   Sum Sq   Mean Sq F value  Pr(>F)  
# matrix(parkoh_prukoh_QTLeffects$est, ncol = 2)[, 2]  1 0.024014 0.0240144  4.9889 0.07582 .
# Coefficients:
# (Intercept)  matrix(parkoh_prukoh_QTLeffects$est, ncol = 2)[, 2]  
# 0.003201                                             1.176388 

# anova(lm(matrix(parkoh_prukoh_QTLeffects$PVE,ncol=2)[,1]~matrix(parkoh_prukoh_QTLeffects$PVE,ncol=2)[,2])) 


# parkon_prukoh_QTLeffects<-rbind(parkon_QTLeffects[c(1,3,4,5,7,8,9),c("est","PVE","QTL")],prukoh_QTLeffects[c(1,3,4,5,6,7,8),c("est","PVE","QTL")])
# anova(lm(matrix(parkon_prukoh_QTLeffects$est,ncol=2)[,2]~matrix(parkon_prukoh_QTLeffects$est,ncol=2)[,1])) 

# #                                                    Df    Sum Sq   Mean Sq F value Pr(>F)  
# #matrix(parkon_prukoh_QTLeffects$est, ncol = 2)[, 1]  1 0.0071240 0.0071240  4.1116 0.0984 .
# #Residuals                                            5 0.0086632 0.0017326                 


# parkon_cereuk_QTLeffects<-rbind(parkon_QTLeffects[c(1,4,5,6,8,9),c("est","PVE","QTL")],cereuk_QTLeffects[c(1,3,4,5,6,7),c("est","PVE","QTL")])
# anova(lm(matrix(parkon_cereuk_QTLeffects$est,ncol=2)[,2]~matrix(parkon_cereuk_QTLeffects$est,ncol=2)[,1])) 

#                                                    Df    Sum Sq   Mean Sq F value Pr(>F)
#matrix(parkon_cereuk_QTLeffects$est, ncol = 2)[, 1]  1 0.0009651 0.0009651     0.2 0.6779
#Residuals                                            4 0.0192995 0.0048249


# prukoh_cereuk_QTLeffects<-rbind(prukoh_QTLeffects[c(1,4,5,7,8),c("est","PVE","QTL")],cereuk_QTLeffects[c(1,3,4,6,7),c("est","PVE","QTL")])
# anova(lm(matrix(prukoh_cereuk_QTLeffects$est,ncol=2)[,2]~matrix(prukoh_cereuk_QTLeffects$est,ncol=2)[,1])) 

# #                                                    Df    Sum Sq   Mean Sq F value   Pr(>F)   
# #matrix(prukoh_cereuk_QTLeffects$est, ncol = 2)[, 1]  1 0.0192304 0.0192304  114.07 0.001755 **
# #Residuals                                            3 0.0005058 0.0001686

# # plot

# cairo_pdf("parkon_prukoh_QTLeffectsize_correlations.pdf", width=4, height=4)
# plot(matrix(parkon_prukoh_QTLeffects$est,ncol=2), xlim=c(0,.25),ylim=c(0,.25), xlab="Effect size (pps) KonPar", ylab="Effect size (pps) PruKoh")
# abline(lm(matrix(parkon_prukoh_QTLeffects$est,ncol=2)[,2]~matrix(parkon_prukoh_QTLeffects$est,ncol=2)[,1]))
# dev.off()

# cairo_pdf("parkon_cereuk_QTLeffectsize_correlations.pdf", width=4, height=4)
# plot(matrix(parkon_cereuk_QTLeffects$est,ncol=2), xlim=c(0,.25),ylim=c(0,.25), xlab="Effect size (pps) KonPar", ylab="Effect size (pps) CerEuk")
# abline(lm(matrix(parkon_cereuk_QTLeffects$est,ncol=2)[,2]~matrix(parkon_cereuk_QTLeffects$est,ncol=2)[,1]))
# dev.off()

# cairo_pdf("prukoh_cereuk_QTLeffectsize_correlations.pdf", width=4, height=4)
# plot(matrix(prukoh_cereuk_QTLeffects$est,ncol=2), xlim=c(0,.25),ylim=c(0,.25), xlab="Effect size (pps) KonPar", ylab="Effect size (pps) CerEuk")
# abline(lm(matrix(prukoh_cereuk_QTLeffects$est,ncol=2)[,2]~matrix(prukoh_cereuk_QTLeffects$est,ncol=2)[,1]))
# dev.off()


# shared_qtl<-data.frame(estimate=c(parkon_QTLeffects$est,prukoh_QTLeffects$est,cereuk_QTLeffects$est),
# QTL = c(parkon_QTLeffects$QTL,prukoh_QTLeffects$QTL,cereuk_QTLeffects$QTL),
# position=c(parkon_QTLeffects$markerpos,prukoh_QTLeffects$markerpos,cereuk_QTLeffects$markerpos),
# cross = c(rep("parkon",nrow(parkon_QTLeffects)),rep("prukoh",nrow(prukoh_QTLeffects)),rep("cereuk",nrow(cereuk_QTLeffects))),
# shared= c("double","unique","triple","triple","triple","double","double","triple","double","triple","unique","triple","triple","triple","double","triple","double","triple","triple","triple","triple","double","triple","double")

# )

# cairo_pdf("Shared_QTL_effects.pdf", width=6.5, height=4)
# ggplot(shared_qtl, aes(x=factor(shared,levels=c("unique","double","triple")),y=estimate,color=as.factor(cross))) + geom_point()
# dev.off()
# aggregate(shared_qtl$estimate, by=list(shared_qtl$shared),FUN=mean)
  # shared          x
# 1  double 0.08512728
# 2  triple 0.08790116
# 3  unique 0.03584714


# TukeyHSD(aov(estimate ~ cross + shared , data=shared_qtl))
  # Tukey multiple comparisons of means
    # 95% family-wise confidence level

# Fit: aov(formula = estimate ~ cross + shared, data = shared_qtl)

# $cross
                      # diff         lwr        upr     p adj
# parkon-cereuk -0.012375189 -0.08736986 0.06261948 0.9081375
# prukoh-cereuk  0.004113967 -0.07290405 0.08113198 0.9899048
# prukoh-parkon  0.016489155 -0.05582093 0.08879925 0.8326785

# $shared
                     # diff        lwr        upr     p adj
# triple-double  0.00277388 -0.0691096 0.07465736 0.9947171
# unique-double -0.04515508 -0.1503817 0.06007156 0.5314976
# unique-triple -0.04792896 -0.1420465 0.04618861 0.4157307

# pheno_dist<-c(1.4,1.6,2.0,3.1)
# avg_effect<-c(0.07, 0.08, 0.09,0.18)

# cairo_pdf("CerEuk_phenodist.pdf",width=4, height=3)
# ggplot(CerEuk_LPmerge$pheno, aes(x=pr)) + geom_histogram(binwidth=0.1, color="black", fill="white") + scale_x_continuous(limits=c(0,4)) + scale_y_continuous(limits=c(0,60))
# dev.off()

# # parkoh_Q1<-find.flanking(parkoh_jitter,5,31) # marker under qtl peak and flanking markers
# # parkon_parkoh_Q1<-find.markerpos(parkon_jitter,marker= markernames(parkon_jitter)[grep(pattern=substr(parkoh_Q1$close,1,7),markernames(parkon_jitter))])
# # prukoh_parkoh_Q1<-find.markerpos(prukoh_jitter,marker= markernames(prukoh_jitter)[grep(pattern=substr(parkoh_Q1$left,1,7),markernames(prukoh_jitter))])

# # parkoh_Q2<-find.flanking(parkoh_jitter,2,70) # marker under qtl peak and flanking markers
# # parkon_parkoh_Q2<-find.markerpos(parkon_jitter,marker= markernames(parkon_jitter)[grep(pattern=substr(parkoh_Q2$right,1,7),markernames(parkon_jitter))])
# # #prukoh_parkoh_Q2<-find.markerpos(prukoh_jitter,marker= markernames(prukoh_jitter)[grep(pattern=substr(parkoh_Q2$close,1,7),markernames(prukoh_jitter))])

# parkoh_Q3<-find.flanking(parkoh_jitter,4,38) # marker under qtl peak and flanking markers
# parkon_parkoh_Q3<-find.markerpos(parkon_jitter,marker= markernames(parkon_jitter)[grep(pattern=substr(parkoh_Q3$left,1,7),markernames(parkon_jitter))])
# prukoh_parkoh_Q3<-find.markerpos(prukoh_jitter,marker= markernames(prukoh_jitter)[grep(pattern=substr(parkoh_Q3$left,1,7),markernames(prukoh_jitter))])

# parkoh_Q4<-find.flanking(parkoh_jitter,1,86) # marker under qtl peak and flanking markers
# parkon_parkoh_Q4<-find.markerpos(parkon_jitter,marker= markernames(parkon_jitter)[grep(pattern=substr(parkoh_Q4$close,1,7),markernames(parkon_jitter))])
# prukoh_parkoh_Q4<-find.markerpos(prukoh_jitter,marker= markernames(prukoh_jitter)[grep(pattern=substr(parkoh_Q4$close,1,7),markernames(prukoh_jitter))])

# parkoh_Q5<-find.flanking(parkoh_jitter,"X",41) # marker under qtl peak and flanking markers
# parkon_parkoh_Q5<-find.markerpos(parkon_jitter,marker= markernames(parkon_jitter)[grep(pattern=substr(parkoh_Q5$right,1,7),markernames(parkon_jitter))])
# prukoh_parkoh_Q5<-find.markerpos(prukoh_jitter,marker= markernames(prukoh_jitter)[grep(pattern=substr(parkoh_Q5$right,1,7),markernames(prukoh_jitter))])

# parkoh_Q7<-find.flanking(parkoh_jitter,6,57) # marker under qtl peak and flanking markers
# parkon_parkoh_Q7<-find.markerpos(parkon_jitter,marker= markernames(parkon_jitter)[grep(pattern=substr(parkoh_Q7$close,1,7),markernames(parkon_jitter))])
# prukoh_parkoh_Q7<-find.markerpos(prukoh_jitter,marker= markernames(prukoh_jitter)[grep(pattern=substr(parkoh_Q7$close,1,7),markernames(prukoh_jitter))])

# parkoh_Q8<-find.flanking(parkoh_jitter,3,74) # marker under qtl peak and flanking markers
# parkon_parkoh_Q8<-find.markerpos(parkon_jitter,marker= markernames(parkon_jitter)[grep(pattern=substr(parkoh_Q8$left,1,7),markernames(parkon_jitter))])
# prukoh_parkoh_Q8<-find.markerpos(prukoh_jitter,marker= markernames(prukoh_jitter)[grep(pattern=substr(parkoh_Q8$left,1,7),markernames(prukoh_jitter))])

# cairo_pdf("Effect_QTL1.pdf", width = 12, height=4)
# par(mfrow=c(1,3))
# plot.pxg(parkon_jitter, marker=parkon_qtl1, pheno.col="pr", ylab="pulse rate", main=paste0("ParKon QTL1 Effect",":",parkon_qtl1), ylim=c(0.5,3.5))
# plot.pxg(parkoh_jitter, marker=parkoh_qtl1, pheno.col="pr", ylab="pulse rate", main=paste0("ParKoh QTL1 Effect",":",parkoh_qtl1), ylim=c(0.5,3.5))
# plot.pxg(prukoh_jitter, marker=prukoh_qtl1, pheno.col="pr", ylab="pulse rate", main=paste0("PruKoh QTL1 Effect",":",prukoh_qtl1), ylim=c(0.5,3.5))
# dev.off()

# cairo_pdf("Effect_QTL2.pdf", width = 12, height=4)
# par(mfrow=c(1,3))
# plot.pxg(parkon_jitter, marker=parkon_qtl2, pheno.col="pr", ylab="pulse rate", main=paste0("ParKon QTL2 Effect",":",parkon_qtl2), ylim=c(0.5,3.5))
# plot.pxg(parkoh_jitter, marker=parkoh_qtl2, pheno.col="pr", ylab="pulse rate", main=paste0("ParKoh QTL2 Effect",":",parkoh_qtl2), ylim=c(0.5,3.5))
# plot.pxg(prukoh_jitter, marker=prukoh_qtl2, pheno.col="pr", ylab="pulse rate", main=paste0("PruKoh QTL2 Effect",":",prukoh_qtl2), ylim=c(0.5,3.5))
# dev.off()

# cairo_pdf("Effect_QTL4.pdf", width = 8, height=4)
# par(mfrow=c(1,2))
# plot.pxg(parkon_jitter, marker=parkon_qtl4, pheno.col="pr", ylab="pulse rate", main=paste0("ParKon QTL4 Effect",":",parkon_qtl4), ylim=c(0.5,3.5))
# plot.pxg(parkoh_jitter, marker=parkoh_qtl4, pheno.col="pr", ylab="pulse rate",  main=paste0("ParKoh QTL4 Effect",":",parkoh_qtl4), ylim=c(0.5,3.5))
# dev.off()

# cairo_pdf("Effect_QTL5b.pdf", width = 12, height=4)
# par(mfrow=c(1,3))
# plot.pxg(parkon_jitter, marker=parkon_qtl5b, pheno.col="pr", ylab="pulse rate", main=paste0("ParKon QTL5b Effect",":",parkon_qtl5b), ylim=c(0.5,3.5))
# plot.pxg(parkoh_jitter, marker=parkoh_qtl5b, pheno.col="pr", ylab="pulse rate",  main=paste0("ParKoh QTL5b Effect",":",parkoh_qtl5b), ylim=c(0.5,3.5))
# plot.pxg(prukoh_jitter, marker=prukoh_qtl5b, pheno.col="pr", ylab="pulse rate",  main=paste0("PruKoh QTL5b Effect",":",prukoh_qtl5b), ylim=c(0.5,3.5))
# dev.off()


# cairo_pdf("Effect_QTL5a.pdf", width = 8, height=4)
# par(mfrow=c(1,2))
# plot.pxg(parkon_jitter, marker=parkon_qtl5a, pheno.col="pr", ylab="pulse rate", main=paste0("ParKon QTL5a Effect",":",parkon_qtl5a), ylim=c(0.5,3.5))
# plot.pxg(parkoh_jitter, marker=parkoh_qtl5a, pheno.col="pr", ylab="pulse rate", main=paste0("ParKoh QTL5a Effect",":",parkoh_qtl5a), ylim=c(0.5,3.5))
# dev.off()

# cairo_pdf("Effect_QTL6.pdf", width = 4, height=4)
# plot.pxg(prukoh_jitter, marker=prukoh_qtl6, pheno.col="pr", ylab="pulse rate", main=paste0("PruKoh QTL6 Effect",":",prukoh_qtl6), ylim=c(0.5,3.5))
# dev.off()

# cairo_pdf("Effect_QTLX.pdf", width = 8, height=4)
# par(mfrow=c(1,2))
# plot.pxg(parkon_jitter, marker=parkon_qtlX, pheno.col="pr", ylab="pulse rate", main=paste0("ParKon QTLX Effect",":",parkon_qtlX), ylim=c(0.5,3.5))
# plot.pxg(parkoh_jitter, marker=parkoh_qtlX, pheno.col="pr", ylab="pulse rate", main=paste0("ParKoh QTLX Effect",":",parkoh_qtlX), ylim=c(0.5,3.5))
# dev.off()



# ## plot QTL plots with Bayes interval:

# #ParKon<-read.delim("QTL_ParKon_HK.txt")
# #ParKoh<-read.delim("QTL_ParKoh_HK.txt")
# #PruKoh<-read.delim("QTL_PruKoh_HK.txt")
# df<-read.delim("combinedQTL_HK_BayesInt.txt")
# #df<-rbind(ParKoh,ParKon,PruKoh)
# #df<-na.omit(df)
# #df$cross<-factor(df$cross)
# df$cross2<-c(rep(2,length(df[df[,"cross"]=="parkoh","cross"])),rep(3,length(df[df[,"cross"]=="prukoh","cross"])),rep(1,length(df[df[,"cross"]=="parkon","cross"])))


# theme_set(theme_bw()+theme(legend.position="none"))


# BayesInt_LK1<-ggplot(df[df[,"chr"]=="1",]) + geom_point(aes(x=cross2, y=pos)) + geom_line(aes(x=cross2,y=pos,group=as.factor(contig)), color="red") + geom_point(aes(x=cross2,y=pos, color=BayesInt))+ geom_segment(aes(x=cross2-0.25, y=QTLcenter, xend=cross2,yend=QTLcenter),arrow = arrow(length = unit(0.03, "npc"))) + scale_y_continuous(limits=c(0,200))
# BayesInt_LK2<-ggplot(df[df[,"chr"]=="2",]) + geom_point(aes(x=cross2, y=pos)) + geom_line(aes(x=cross2,y=pos,group=as.factor(contig)), color="red") + geom_point(aes(x=cross2,y=pos, color=BayesInt))+ geom_segment(aes(x=cross2-0.25, y=QTLcenter, xend=cross2,yend=QTLcenter),arrow = arrow(length = unit(0.03, "npc"))) + scale_y_continuous(limits=c(0,260))
# BayesInt_LK3<-ggplot(df[df[,"chr"]=="3",]) + geom_point(aes(x=cross2, y=pos)) + geom_line(aes(x=cross2,y=pos,group=as.factor(contig)), color="red") + geom_point(aes(x=cross2,y=pos, color=BayesInt))+ geom_segment(aes(x=cross2-0.25, y=QTLcenter, xend=cross2,yend=QTLcenter),arrow = arrow(length = unit(0.03, "npc"))) + scale_y_continuous(limits=c(0,200))
# BayesInt_LK4<-ggplot(df[df[,"chr"]=="4",]) + geom_point(aes(x=cross2, y=pos)) + geom_line(aes(x=cross2,y=pos,group=as.factor(contig)), color="red") + geom_point(aes(x=cross2,y=pos, color=BayesInt))+ geom_segment(aes(x=cross2-0.25, y=QTLcenter, xend=cross2,yend=QTLcenter),arrow = arrow(length = unit(0.03, "npc"))) + scale_y_continuous(limits=c(0,200))
# BayesInt_LK5<-ggplot(df[df[,"chr"]=="5",]) + geom_point(aes(x=cross2, y=pos)) + geom_line(aes(x=cross2,y=pos,group=as.factor(contig)), color="red") + geom_point(aes(x=cross2,y=pos, color=BayesInt))+ geom_segment(aes(x=cross2-0.25, y=QTLcenter, xend=cross2,yend=QTLcenter),arrow = arrow(length = unit(0.03, "npc"))) + scale_y_continuous(limits=c(0,200))
# BayesInt_LK6<-ggplot(df[df[,"chr"]=="6",]) + geom_point(aes(x=cross2, y=pos)) + geom_line(aes(x=cross2,y=pos,group=as.factor(contig)), color="red") + geom_point(aes(x=cross2,y=pos, color=BayesInt))+ geom_segment(aes(x=cross2-0.25, y=QTLcenter, xend=cross2,yend=QTLcenter),arrow = arrow(length = unit(0.03, "npc"))) + scale_y_continuous(limits=c(0,200))
# BayesInt_LK7<-ggplot(df[df[,"chr"]=="7",]) + geom_point(aes(x=cross2, y=pos)) + geom_line(aes(x=cross2,y=pos,group=as.factor(contig)), color="red") + geom_point(aes(x=cross2,y=pos, color=BayesInt))+ geom_segment(aes(x=cross2-0.25, y=QTLcenter, xend=cross2,yend=QTLcenter),arrow = arrow(length = unit(0.03, "npc"))) + scale_y_continuous(limits=c(0,200))
# BayesInt_LKX<-ggplot(df[df[,"chr"]=="X",]) + geom_point(aes(x=cross2, y=pos)) + geom_line(aes(x=cross2,y=pos,group=as.factor(contig)), color="red") + geom_point(aes(x=cross2,y=pos, color=BayesInt))+ geom_segment(aes(x=cross2-0.25, y=QTLcenter, xend=cross2,yend=QTLcenter),arrow = arrow(length = unit(0.03, "npc"))) + scale_y_continuous(limits=c(0,200))

# cairo_pdf("combinedQTL_BayesInt.pdf",width=4,height=32)
# grid.arrange(BayesInt_LK1,BayesInt_LK2,BayesInt_LK3,BayesInt_LK4,BayesInt_LK5,BayesInt_LK6,BayesInt_LK7,BayesInt_LKX,ncol=1)
# dev.off()


