#THIS CODE IS FOR THE GOBY ADULT SEX RATIO MANUSCRIPT
#Extracting the list of sex-associated markers from the outputs of Stacks and RADsex


#setwd("C://Users//ivainm//AppData//Local//Packages//CanonicalGroupLimited.Ubuntu18.04LTS_79rhkp1fndgsc//LocalState//rootfs//home//ivainm//realdir//Radsex_files")
setwd("C://Users//ivmar9032//Documents//DYNAMAR local//Manuscripts//sex_det//data_analysis")


#Packages----
#install.packages("devtools")
#devtools::install_github("SexGenomicsToolkit/sgtr")

library(sgtr)
library(ggplot2)
library(ggExtra)

#_____-----
# RADsex markers  ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Haplotype 1 of the reference genome
#radmaps<-read.table("mapsex_hap1.tsv", header=T)
#radhits<-radmaps[radmaps$Signif=="True",]
#write_tsv(radhits,"radsex_hap1.tsv")
radhits<-read.table("radsex_hap1.tsv", header=T)
str(radhits)

#Haplotype 2 of the reference genome
#radmaps2<-read.table("mapsex_hap2.tsv", header=T)
#radhits2<-radmaps2[radmaps2$Signif=="True",]
#write_tsv(radhits2,"radsex_hap2.tsv")
radhits2<-read.table("radsex_hap2.tsv", header=T)

### Plot radsex markers----
radsex_distrib("radsex_hap2.tsv", groups = c("m", "f"), group_labels = c("Males", "Females"), output_file = "figures_forward/radsex_distrib1_1.png")
radsex_distrib("distribution1_3.tsv", groups = c("m", "f"), group_labels = c("Males", "Females"), output_file = "figures_forward/radsex_distrib1_3.png")
radsex_distrib("distribution1_5.tsv", groups = c("m", "f"), group_labels = c("Males", "Females"), output_file = "figures_forward/radsex_distrib1_5.png")
radsex_distrib("distribution1_10.tsv", groups = c("m", "f"), group_labels = c("Males", "Females"), output_file = "figures_forward/radsex_distrib1_10.png")


# SNPs from Stacks, allele frequency is associated with sex ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Haplotype 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assochap1<-read.table("assocsex_hap1.assoc.fisher", header =T)
str(assochap1)
subass<-assochap1[assochap1$P<0.001,] #avoid plotting everything, subset by p value
str(subass)

ggplot(data=subass[subass$P<1,], aes(x=F_A, y=F_U))+ #find a satisfying p value threshold using the plot
  geom_point( aes(fill=P),shape=21, size=3, alpha=0.5)+
  labs(x="Minor allele frequency in males",y="Minor allele frequency in females",title="Association of allele frequency with sex")+
  scale_fill_gradientn(colours = c("red", "orange", "lightblue", "royalblue"),
                       values = c(0,0.0000000001, 0.001, 0.01,  1))

#Haplotype 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assochap2<-read.table("assocsex_hap2.assoc.fisher", header =T)
str(assochap2)
subass2<-assochap2[assochap2$P<1,]
str(subass2)

ggplot(data=subass2[subass2$P<0.0000000001,], aes(x=F_A, y=F_U))+ #use same p value threshold as for hap1
  geom_point( aes(fill=P),shape=21, size=3)+
  labs(x="Minor allele frequency in males",y="Minor allele frequency in females",title="Association of allele frequency with sex")+
  scale_fill_gradientn(colours = c("darkblue", "blue", "lightblue", "red"),
                       values = c(0,0.001, 0.05, 0.1,  1))


#Using the p value threshold defined graphically above, export list of markers
candfreq_hap1_m<-subass[subass$P<0.0000000001 &  #markers polymorphic in males, haplotype 1
                          subass$F_A>0.3,]
dim(candfreq_hap1_m)
candfreq_hap1_f<-subass[subass$P<0.0000000001 &  #markers polymorphic in females, hap1
                          subass$F_U>0.5,]
candfreq_hap2_m<-subass2[subass2$P<0.0000000001 & #markers polymorphic in males, haplotype 2
                           subass2$F_A>0.3,]
dim(candfreq_hap1_m)
str(candfreq_hap2_m)
candfreq_hap2_f<-subass2[subass2$P<0.0000000001 & #markers polymorphic in females, hap2
                           subass2$F_U>0.5,]
length(unique(rbind(candfreq_hap1_m,candfreq_hap2_m)$SNP))

## export marker list----
write_tsv(candfreq_hap1_m,"stacks_freq_hap1_m.txt")
write_tsv(candfreq_hap2_m,"stacks_freq_hap2_m.txt")


# SNPs from Stacks, missingness is associated with sex ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Haplotype 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
missinghap1<-read.table("missingsex_hap1.missing", header =T)
submiss<-missinghap1[missinghap1$P<0.05,]
str(submiss)

ggplot(data=submiss[submiss$P<0.000001,], aes(x=F_MISS_A, y=F_MISS_U))+ #define p value threshold graphically
  geom_point( aes(fill=P),shape=21, size=3)+
  labs(x="Missing call frequency in males",y="Missing call frequency in females",title="Association of locus missingness with sex")+
  scale_fill_gradientn(colours = c("darkblue", "blue", "lightblue", "red"),
                       values = c(0,0.001, 0.05, 0.1,  1))
#Haplotype 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
missinghap2<-read.table("missingsex_hap2.missing", header =T)
str(missinghap2)
submiss2<-missinghap2[missinghap2$P<0.05,]
str(submiss2)


ggplot(data=submiss2[submiss2$P<0.000001,], aes(x=F_MISS_A, y=F_MISS_U))+ #use same p value threshold 
  geom_point( aes(fill=P),shape=21, size=3)+
  labs(x="Missing call frequency in males",y="Missing call frequency in females",title="Association of locus missingness with sex")+
  scale_fill_gradientn(colours = c("darkblue", "blue", "lightblue", "red"),
                       values = c(0,0.001, 0.05, 0.1,  1))



#Subset data to extract position of SNPs of interest and write file
candmiss_hap1_m<-submiss[submiss$P<0.000001 &
                           submiss$F_MISS_U>0.25,]
candmiss_hap1_f<-submiss[submiss$P<0.000001 &
                           submiss$F_MISS_A>0.4,]
candmiss_hap2_m<-submiss2[submiss2$P<0.000001 &
                            submiss2$F_MISS_U>0.25,]
candmiss_hap2_f<-submiss2[submiss2$P<0.000001 &
                            submiss2$F_MISS_A>0.4,]

## export marker list----
write_tsv(candmiss_hap1_m,"stacks_miss_hap1_m.txt")
write_tsv(candmiss_hap1_f,"stacks_miss_hap1_f.txt")
write_tsv(candmiss_hap2_m,"stacks_miss_hap2_m.txt")
write_tsv(candmiss_hap2_f,"stacks_miss_hap2_f.txt")



#____----
#  Genotype of Adult males and females at the sex markers----
setwd("C://Users//ivmar9032//Documents//DYNAMAR local//Manuscripts//sex_det//data_analysis//sexID")
library(vcfR)

#vcf file from ddRAD, aligned to minutus
vcf <- read.vcfR("vcf_files/populations.snps.vcf", verbose = FALSE)
#rm(vcf)
males<-read.table("vcf_files/maleIDs.txt", header=F) # Corrected for plate swap 7&8.
females<-read.table("vcf_files/femaleIDs.txt", header=F) 
#List of sex markers from ddRAD dataset
sexmarkers <- read.table("Sex_loci_RAD.txt", header = TRUE) %>% mutate(Site = paste0(CHR, "_", "POS", "_"))
## Filter for the sex markers
# Extract the chromosome and position data from the VCF file
vcf_chrom <- getCHROM(vcf) # Chromosome information
vcf_pos <- getPOS(vcf)     # Position information

# Merge chromosome and position information into a data frame
vcf_sites <- data.frame(CHR = vcf_chrom, POS = vcf_pos)

# Create a logical vector to match rows in vcf_sites with those in sexmarkers
# This will be TRUE for rows in vcf that match chromosome and position in sexmarkers
matching_sites <- with(vcf_sites, (CHR %in% sexmarkers$CHR) & (POS %in% sexmarkers$POS))

# Now subset the VCF file to include only the matching sites
vcfsex <- vcf[matching_sites, ]

# Get genotypes (I haven't used this so far)
gt <- extract.gt(vcfsex, element='GT', return.alleles = TRUE) %>% data.frame()
gt <- tibble::rownames_to_column(gt, "Site") %>% data.frame()
gt <- reshape2::melt(gt, id = c("Site")) 
names(gt) <- c("Site", "Sample", "alleles")
gt <- gt %>% tidyr::separate(alleles, c("Genotype1", "Genotype2"), sep = "/", remove = TRUE)

## Get the allele depth
ad <- extract.gt(vcfsex, element='AD') %>% data.frame()
ad <- tibble::rownames_to_column(ad, "Site") %>% data.frame()

# Reshape into long format
ad <- reshape2::melt(ad, id = c("Site"))
names(ad) <- c("Site", "Sample", "AD")
# Get the different alelle depts as separate columns
ad <- ad %>% tidyr::separate(AD, c("Allele1", "Allele2", "Allele3"), sep = ",", remove = TRUE)
ad$Allele1 <- as.numeric(ad$Allele1)
ad$Allele2 <- as.numeric(ad$Allele2)
ad$Allele3 <- as.numeric(ad$Allele3)
# Often the there is no third allele, so fix that to be 0 instead of NA
ad <- ad %>% mutate(Allele3 = case_when(is.na(Allele3) ~ 0, TRUE ~ Allele3)) # Fix the third allele

# Calculate their allele frequencies
ad <- ad %>% rowwise() %>% mutate(DP = Allele1 + Allele2 + Allele3) %>% mutate(freq1 = Allele1/DP, freq2 = Allele2/DP, freq3 = Allele3/DP) 

# Calculate the major (big) allele frequency 
MAFs <- ad %>% rowwise() %>% mutate(baf = max(Allele1, Allele2, Allele3)/DP) %>% data.frame() # rowwise makes the trick to do it per row
# Add the genotypes (I didn't use this in the end)
MAFs <- merge(MAFs, gt, by = c("Site", "Sample"))

## Add metadata
MAFs$SEX <-NA
MAFs[MAFs$Sample %in% males$V1,]$SEX<-"male"
MAFs[MAFs$Sample %in% females$V1,]$SEX<-"female"
table(MAFs$SEX)
table(MAFs$Site)

## identify Hm and Hz individuals
ggplot(MAFs,aes(x=freq1))+geom_histogram()+
  facet_grid(rows=vars(SEX),cols=vars(Site))

MAFs$Heterozygote<-NA
MAFs$Homozygote<-NA
MAFsub<-MAFs[which(!is.na(MAFs$freq1)),]
MAFsub[MAFsub$freq1>0.9,]$Homozygote<-"yes"
MAFsub[MAFsub$freq1<0.75 & MAFsub$freq1>0.25,]$Heterozygote<-"yes"


# format summary tables
Hzs<-data.frame(table(MAFsub$SEX, MAFsub$Heterozygote,MAFsub$Site))
Hms<-data.frame(table(MAFsub$SEX, MAFsub$Homozygote,MAFsub$Site))


Hz<-cbind(Hzs[Hzs$Var1=="female",],Male=Hzs[Hzs$Var1=="male","Freq"])
Hz<-Hz[,3:5]
colnames(Hz)<-c("site","Hz.females","Hz.males")
Hz

Hm<-cbind(Hms[Hms$Var1=="female",],Male=Hms[Hms$Var1=="male","Freq"])
Hm<-Hm[,3:5]
colnames(Hm)<-c("site","Hm.females","Hm.males")
Hm

#calculate proportions of heterozygotes and homozygotes males and females
glob<-cbind(Hm,Hz[,2:3])
glob$mHZperc<-glob$Hz.males/(glob$Hz.males+glob$Hm.males)
glob$fHMperc<-glob$Hm.females/(glob$Hm.females+glob$Hz.females)
glob
mean(glob$mHZperc)
mean(glob$fHMperc)