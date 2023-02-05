#UVMR

exp_dat <- read_exposure_data(
  filename = 'data1.csv',
  clump = FALSE,
  sep= "",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col ="Tested_Allele",
  other_allele_col = "Other_Allele",
  eaf_col = "Freq_Tested_Allele",
  pval_col = "P",  pos_col = "pos",chr_col = "chr",ncase_col = "N",
)
install.packages("readxl")
library(readxl)
data1 <- read_excel( 'Waist-to-hip ratio.xlsx',sheet = 1) 

#online 
exp_dat <- extract_instruments(
  outcomes='ukb-b-9405',access_token= NULL
)
dim(exp_dat)


#local
setwd("C:/Users/Administrator/Desktop/pone")
data1=read.table("ieu-b-40.vcf",sep="\t",header = T)
data2=read.table("WC.csv",sep=",",header = T)

data2=read.table("NamjouB_31311600_NAFLD.txt",sep="\t",header = T)
data1=read.xlsx("Mets.xlsx",sep="",header = T,sheetIndex = 1)
data1$Phenotype <- 'nafld'   
exp_dat=format_data(
  data2,
  snps = NULL, type='exposure',
  header = TRUE,
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value"
)

dim(exp_dat)
exp_dat=data2
exposure_dat=subset(exp_dat, exp_dat$pval.exposure
                    <5e-08)
dim(exposure_dat)

exp_dat <-clump_data(exposure_dat,clump_r2=,clump_kb=)


#local outcome
dim(exp_dat)
data2=read.csv("liver.csv",sep=",",header = T)
data1=read.table("Ghodsian_meta_NAFLD.tsv",sep="\t",header = T)
data2$Phenotype <- 'liver' 
out_dat=format_data(
  data1,
  snps = NULL, type='outcome',
  header = TRUE,
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value"
)


#online 



bc <- extract_outcome_data(
  snps=exp_dat$SNP,
  outcomes='',
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL
)
dim(bc)

mydata <- harmonise_data(
  exposure_dat=exp_dat,
  outcome_dat=bc,
  action= 2
)

mr_rucker_cooksdistance(mydata, parameters = default_parameters())


res <- mr(mydata)
res
het <- mr_heterogeneity(mydata)
het

mr(mydata,method_list=c('mr_ivw_mre'))


pleio <- mr_pleiotropy_test(mydata)
pleio

single <- mr_leaveoneout(mydata)
p1=mr_leaveoneout_plot(single)

p2=mr_scatter_plot(res,mydata)


res_single <- mr_singlesnp(mydata)
p3=mr_forest_plot(res_single)

p4=mr_funnel_plot(res_single)












#MVMR

library(forestmodel)
library(TwoSampleMR)
id_exposure <- c("","") 
id_outcome <- ""
data3=read.table("",sep=",",header = T)
exposure_dat=data3#

exposure_dat <- mv_extract_exposures(id_exposure)

mvmr=read.csv("",sep=",",header = T) #local

mvmr=clump_data(exposure_dat,clump_r2 = 0.01,clump_kb = 1000)
dim(mvmr)

outcome_dat <- extract_outcome_data(mvmr$SNP, id_outcome) 
mvdat <- mv_harmonise_data(mvmr,  outcome_dat) # harmonisation
res <- mv_multiple(mvdat) #1
res=mv_residual(mvdat)#2
res
