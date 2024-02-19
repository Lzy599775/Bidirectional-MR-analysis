R code for Bidirectional-MR-analysis

exposureName=" "         
outcomeName=" "          
inputFile=" "            
outcomeFile=" "         
samplesize.exposure<- 
samplesize.outcome <-
    

tsv1<- vroom::vroom(inputFile)
names(tsv1)
exposuredata<-TwoSampleMR::format_data(tsv1,type = "exposure",
                                       phenotype_col = "Phenotype",
                                       snp_col = "rsid",beta_col = "BETA",
                                       se_col = "SE",effect_allele_col = "Allele1",
                                       other_allele_col = "Allele2",pval_col = "p_value",
                                       eaf_col = "Freq1",chr_col = "chromosome",
                                       pos_col = "base_pair_location",
                                       samplesize_col = "Samplesize.exposure")
exp.P<-subset(exposuredata, pval.exposure<5e-08)
write.csv(exp.P, file="exp.p.csv", row.names=T)

tsv2<- vroom::vroom(outcomeFile)
names(tsv2)
outcomeData<-TwoSampleMR::format_data(tsv2,type = "outcome",
                                      phenotype_col = "Phenotype",
                                      snp_col = "rsids",beta_col = "beta",
                                      se_col = "sebeta",effect_allele_col = "alt",
                                      other_allele_col = "ref",eaf_col = "af_alt",
                                      samplesize_col = "samplesize.outcome",
                                      pval_col = "pval",chr_col = "#chrom",pos_col = "pos")
write.csv(outcomeData, file="out.csv", row.names=T)

exposureFile="exp.p.csv"    
exposure_dat<-read_exposure_data(filename=exposureFile,
                                 phenotype_col = "exposure",
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 pval_col = "pval.exposure",
                                 se_col = "se.exposure",pos_col = "pos.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 eaf_col = "eaf.exposure",chr_col = "chr.exposure",
                                 samplesize_col = "samplesize.exposure",
                                 clump = F)

exposure_dat_clumped <- clump_data(exposure_dat,
                                   clump_kb=10000, clump_r2=0.001)
write.csv(exposure_dat_clumped, file="exp.clumped.csv", row.names=F)



########## Find and remove confounding factors
inputFile="exp.clumped.csv"
dat=read.csv(inputFile, header=T, sep=",", check.names=F)
snpId=dat$SNP
y=seq_along(snpId)
chunks <- split(snpId, ceiling(y/100))
outTab=data.frame()
outTab <- data.frame()
for (i in names(chunks)) {
    confounder <- phenoscanner(
        snpquery = chunks[[i]],
        catalogue = "GWAS",
        pvalue = 5e-08,
        proxies = "None",
        r2 = 0.8,
        build = 37
    )
    outTab <- rbind(outTab, confounder$results)
}
write.csv(outTab, "confounder.result.csv", row.names=F)

res_snp_confounder2<-outTab
confounder_factor<-c("smok","body mass","educa","obes",
                     "alco","inco","stre","depres")   
confounder_result<-data.frame()

for (factor in confounder_factor) {
    factor_confounder<-res_snp_confounder2[grepl(factor,res_snp_confounder2$trait),]
    confounder_result<-rbind(confounder_result,factor_confounder)
    confounder_result <- rbind(confounder_result, factor_confounder)
}

confounder_snp<-unique(confounder_result$snp)
confounder_snp
de1_SNP<-confounder_snp
delSnp <- c()
dat=dat[!dat$SNP %in% delSnp,]
write.csv(dat, file="exp.confounder.csv", row.names=F)

inputFile<-"exp.confounder.csv"
dat<-read.csv(inputFile, header=T, sep=",", check.names=F)
datR = transform(dat,
                 R2 = 2 * ((beta.exposure)^2 * eaf.exposure * (1 - eaf.exposure)) /
                     (2 * (beta.exposure)^2 * eaf.exposure * (1 - eaf.exposure) +
                          2 * (se.exposure)^2 * samplesize.exposure * eaf.exposure * (1 - eaf.exposure)))

datFR=transform(datR,F=(samplesize.exposure-2)*R2/(1-R2))      
result <- sum(datFR$R2)
Ffilter=10
exp.FR=datFR[datFR$F>Ffilter,]
write.csv(exp.FR, "exp.FR.csv", row.names=F)

exposureFile="exp.FR.csv"
exposure_dat<-read_exposure_data(filename=exposureFile,
                                 phenotype_col = "exposure",
                                 sep = ",",samplesize_col = "samplesize.exposure",
                                 snp_col = "SNP",
                                 pval_col = "pval.exposure",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 eaf_col = "eaf.exposure",
                                 chr_col = "chr.exposure",
                                 pos_col = "pos.exposure",id_col = "id.exposure",
                                 clump = F)
outcomeTab<-merge(exposure_dat, outcomeData, by.x="SNP", by.y="SNP")
write.csv(outcomeTab, file="outcomedat.csv")

outcome_dat<-read_outcome_data(snps=exposure_dat$SNP,
                               filename="outcomedat.csv",
                               phenotype_col = "outcome",
                               sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta.outcome",
                               se_col = "se.outcome",
                               eaf_col = "eaf.outcome",
                               effect_allele_col = "effect_allele.outcome",
                               other_allele_col = "other_allele.outcome",
                               pval_col = "pval.outcome",
                               chr_col = "chr.outcome",
                               pos_col = "pos.outcome",id_col = "id.outcome")
exposure_dat$exposure=exposureName
outcome_dat$outcome=outcomeName
dat2<-harmonise_data(exposure_dat=exposure_dat,outcome_dat=outcome_dat)
dat2$samplesize.outcome<-samplesize.outcome
dat2=dat2[dat2$mr_keep=="TRUE",]
write.csv(dat2, file="SNP.csv", row.names=F)


Steiger_test <- directionality_test(dat2)
Steiger_test
write.csv(Steiger_test,file="Steiger.csv")


mrResult=mr(dat2)   
OR=generate_odds_ratios(mrResult)
write.csv(OR, file="MR-OR.csv", row.names=F)


heterTab=mr_heterogeneity(dat2)
write.csv(heterTab, file="heterogeneity.csv", row.names=F)


mr_pleiotropy_test(dat2) 
pleioTab=mr_pleiotropy_test(dat2)
write.csv(pleioTab, file="pleiotropy.csv", row.names=F)


set.seed(1234)
presso=run_mr_presso(dat2,NbDistribution = 3000)  
presso
rows_to_remove<-c() 
dat2 <- dat2[-rows_to_remove,]
write.csv(dat2, file = "MR-PRESSO.csv")


png(filename = 'funnel.tif',width = 4500,height = 2800,res=350)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat2))
dev.off()

png(filename = 'scatter.tif',width = 4500,height = 2800,res=350)
mr_scatter_plot(mrResult, dat2)
dev.off()

res_single=mr_singlesnp(dat2)   
png(filename = 'forest.tif',width = 4500,height = 2800,res=350)
mr_forest_plot(res_single)
dev.off()

mr_leaveoneout(dat2)
png(file="leaveoneout.tif", width = 4500,height = 2800,res=350)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat2))
dev.off()

res_single <- mr_singlesnp(dat2)
singleSNP=res_single
write.csv(singleSNP, file="singleSNP.csv", row.names=F)

dat2R2=transform(dat2,
                 R2 = 2 * ((beta.exposure)^2 * eaf.exposure * (1 - eaf.exposure)) /
                     (2 * (beta.exposure)^2 * eaf.exposure * (1 - eaf.exposure) +
                          2 * (se.exposure)^2 * 317717 * eaf.exposure * (1 - eaf.exposure)))
result <- sum(dat2R2$R2)
print(result)
