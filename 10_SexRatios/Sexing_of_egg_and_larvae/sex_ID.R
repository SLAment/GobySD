#THIS CODE IS FOR THE GOBY ADULT SEX RATIO MANUSCRIPT
#setwd("C://Users//ivainm//Working_Document//DYNAMAR//popGen//data_analysis//sexID")
setwd("C://Users//ivmar9032//Documents//DYNAMAR local//Manuscripts//sex_det//data_analysis//sexID")

#Packages----
library(vcfR)
library(effects)
library(car)
library(ggplot2)
library(dplyr) 
library(tidyverse)
library(stringr) # for str_remove
library(reshape2) # To fix dataframe with melt

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# _________ ----
#Genotypic sex of eggs and larvae----
#-Identifying genotypic sex of eggs, larvae and unsexed adults
#-Analysis of the sex ratio data from these samples
#-Associated plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##Loading files------------------------------------------------------------------
#Info on sample collection date, area, type (larvae, egg, ref...etc)
popmap<-read.table("popmap_corrected.txt", header=T) # Corrected for plate swap 7&8.
str(popmap)
#how many nests per location are used?
table(popmap$nest, popmap$Area)
#5 nests per location, 35 samples per nest

#Info of %mapping of the reads for each sample
maps<-read.delim("general_stats_table.tsv", header=T, sep="\t") 
maps2<-maps[,c(1,5)] #get rid of useless columns
colnames(maps2)<-c("ID","Percent_Mapped")
maps2$ID<-gsub("SORT-", "", maps2$ID) #format the ID to match other datasets
str(maps2)
#Info on size of the fasta files 
sizes<-read.table("filesize.txt", header=T)
colnames(sizes)<-c("filesize","ID")
sizes$ID<-gsub("_L001_R1_001.fastq.gz", "", sizes$ID)#clean up file ID
sizes$ID<-gsub("_L001_R2_001.fastq.gz", "", sizes$ID)
sizes$ID<-gsub("_S", "", sizes$ID)
sizes<-sizes[!duplicated(sizes$ID),] #each individual has two files because of forward and reverse reads. Keep one.
str(sizes) #sizes are in Kb
#vcf file from GTseq
vcf <- read.vcfR("sexcalls2.vcf", verbose = FALSE)
#List of sex markers from ddRAD dataset
sexmarkers <- read.table("Sex_loci_RAD.txt", header = TRUE) %>% mutate(Site = paste0(CHR, "_", "POS", "_"))


##Filtering of individuals-------------------------------------------------------
# We have decided to filter out samples for which the file size is too small, indicating few reads (<1500K fasta)
ggplot(sizes[sizes$filesize<10000,],aes(x=filesize))+geom_histogram() #the subsetting is to ignore the giant file with all the reads not demultiplexed
#This histogram shows that the cutoff at 1.5Mb is straightforward.
#How many files are below the threshold?
length(sizes[sizes$filesize<1500,1])
#41. let's keep there ID in mind
IDshort<-sizes[sizes$filesize<1500,]$ID

# We also filter out individuals with alignment score <60%
ggplot(maps2, aes(x=Percent_Mapped))+geom_histogram() # 60% seems reasonable. A bit lenient but the size filtering takes care of the doubtful ones anyway
#who are the goners?
length(maps2[maps2$Percent_Mapped<60,1]) #23, but a few are already gone because of size filtering
IDbadmap<-maps2[maps2$Percent_Mapped<60,1] #keep them in mind
maps2[maps2$Percent_Mapped<60 & !(maps2$ID %in% IDshort),1] #who are the new guys that were not already filtered out by size?



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Process vcf (Lore's code) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

# existing_individuals <- factor(colnames(vcf@gt)[-1]) %>% # Exclude the first column which is FORMAT
#   gsub("-", ".", .) # Replace the - for . because the processing later changes the names


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate the MAF (Lore's code) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## This is a bcftools vcf

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
str(MAFs)
## Add metadata
IDs <- str_remove(MAFs$Sample, "SORT.") %>% str_remove(., ".bam")
MAFs$ID <- as.numeric(IDs)
#MAFs <- merge(MAFs, sampleinfo, by = "ID")
MAFs <- merge(MAFs, popmap, by = "ID")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Finding patterns (Lore's code with Ivain's edits) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## What samples have no coverage for any sex-linked marker?
# Make a new data frame for each sample to decide their genetic sex
SampleDF <- MAFs %>% 
  filter(!is.na(baf)) %>%  # remove missing data
  group_by(Sample) %>% summarise(totalDP = sum(DP), 
                                 avgBaf = mean(baf), # average major allele freq
                                 n = n()) %>% # count how many markers made it to the average
  mutate(genetic_sex = case_when(avgBaf > 0.9 ~ "Female", TRUE ~ "Male"))

# Histogram of the total allele depth
ggplot(SampleDF, aes(x = totalDP)) + geom_histogram() 
#We can easily make a cutoff again below 1000 depth. Let's see if we filter out new individuals
length(SampleDF[SampleDF$totalDP<1000,]$Sample)
IDshallow<-gsub("SORT.", "", SampleDF[SampleDF$totalDP<1000,]$Sample)#clean up file ID
IDshallow<-gsub(".bam","",IDshallow)
#who are the new samples that haven't been already filtered out by size or % mapped?
maps2[!(maps2$ID %in% IDbadmap) & !(maps2$ID %in% IDshort) & (maps2$ID %in% IDshallow),1] 
#23 more removed
length(maps2[(maps2$ID %in% IDbadmap) | (maps2$ID %in% IDshort) | (maps2$ID %in% IDshallow),1] ) #total 67 individuals filtered out

# Global average major allele freq distribution
ggplot(SampleDF, aes(x = avgBaf)) + geom_histogram() +
  theme_bw() + xlab("Major allele frequency weighted by number of markers")
# There are two very clear groups

# How does the major allele frequency look like if I divide them with genetic_sex
str(SampleDF) 
#select(SampleDF, Sample, genetic_sex)
#SampleDF[,c("Sample","genetic_sex")]
MAFs <- merge(MAFs, SampleDF[,c("Sample","genetic_sex")], by = "Sample")
#MAFs <- merge(MAFs, select(SampleDF, Sample, genetic_sex), by = "Sample")

### Major allele distribution of the markers for all samples
(majorallplot <- ggplot(MAFs, aes(baf, fill = genetic_sex)) + 
    geom_histogram() +
    facet_grid(genetic_sex~Site) + 
    xlab("Major allele frequency") +
    theme_bw() + theme(legend.position = "none"))

ggsave("majorallplot.pdf", plot = majorallplot, width = 11, height = 3.5)

# scaffold2702_47682 has a lot of missing data
# scaffold7908_69411 and scaffold709_69413 look strange, why is there a bump around 0.75? a paralog?
# scaffold4193_6290 is not exactly pretty but it goes more in the right direction

ggplot(MAFs %>% filter(type %in% c("male_control", "female_control")), aes(baf, fill = genetic_sex)) + 
  geom_histogram() +
  facet_grid(type~Site) + 
  xlab("Major allele frequency") +
  theme_bw() + theme(legend.position = "none")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## PCA of the samples with just the nice markers ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badSNPs <- c("scaffold2702_47682", "scaffold709_69411", "scaffold709_69413", "scaffold4193_6290")

# Re-calculate avgBaf
SampleDFpretty <- MAFs %>% 
  filter(!is.na(baf)) %>% 
  filter(!Site %in% badSNPs) %>% 
  group_by(Sample) %>% summarise(totalDP = sum(DP), 
                                 avgBaf = mean(baf), # average major allele freq
                                 n = n()) %>% # count how many markers made it to the average
  mutate(genetic_sex = case_when(avgBaf > 0.9 ~ "Female", TRUE ~ "Male"))
str(MAFs)
# Make a matrix of the allele frequencies
allele_freqs <- MAFs %>% 
  dplyr::select(Site, Sample, baf) %>%
  dplyr::filter(!Site %in% badSNPs) %>% 
  group_by(Site, Sample) %>%
  filter(!is.na(baf)) %>% # Remove rows with missing data
  tidyr::pivot_wider(names_from = Site, values_from = baf) # Spread the Sites into columns
allele_freqs$SampleID<-gsub("SORT.","",allele_freqs$Sample)
allele_freqs$SampleID<-gsub(".bam","",allele_freqs$SampleID)

# Check and remove any non-numeric columns just in case
allele_data <- allele_freqs[, sapply(allele_freqs, is.numeric)]
allele_data$Sample <- allele_freqs$SampleID # put sample name back


# Remove samples and sites with missing data
allele_data_clean <- allele_data[complete.cases(allele_data), ]

#Alternatively: do a second PCA removing all samples that should be filtered out
allele_data_clean_filtered<-allele_data_clean[!(allele_data_clean$Sample %in% IDshort) &
                                                !(allele_data_clean$Sample %in% IDshallow) &
                                                !(allele_data_clean$Sample %in% IDbadmap) ,]


# Perform PCA #1 (clean, unfiltered)
pca_result <- prcomp(allele_data_clean %>% dplyr::select(!Sample), center = TRUE, scale. = TRUE) 
setdiff(MAFs %>% dplyr::select(Site, Sample, baf) %>% .$Sample, allele_data_clean$Sample) %>% length() # 23 samples were removed (15 have no cov in any marker)

# Extract the PCA scores for the first two components
pca_data <- as.data.frame(pca_result$x)
pca_data %>% dim() # 730 out of 768 samples

# Add sample IDs for labeling
pca_data$Sample <- allele_data_clean$Sample
pca_data$nums <- rownames(allele_data_clean)
# Add more information about the samples
pca_data <- merge(pca_data, SampleDF, by = "Sample")
pca_data <- pca_data %>% mutate(genetic_sex = case_when(avgBaf > 0.9 ~ "Female", TRUE ~ "Male"))

## Plot the PCA 
(PCAp <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = totalDP, shape = genetic_sex)) +
    theme_bw() + geom_point(alpha = 0.7, size = 3) +
    scale_colour_continuous(type = "viridis") )

ggsave("sexmarkersPCA.pdf", plot = PCAp, width = 6, height = 5)

# Perform PCA #2 (clean, and filtered)
pca_result2 <- prcomp(allele_data_clean_filtered %>% dplyr::select(!Sample), center = TRUE, scale. = TRUE) 
setdiff(MAFs %>% dplyr::select(Site, Sample, baf) %>% .$Sample, allele_data_clean_filtered$Sample) %>% length() 

# Extract the PCA scores for the first two components
pca_data2 <- as.data.frame(pca_result2$x)
pca_data2 %>% dim() # 701 out of 768 samples

# Add sample IDs for labeling
pca_data2$Sample <- allele_data_clean_filtered$Sample
pca_data2$nums <- rownames(allele_data_clean_filtered)
# Add more information about the samples
SampleDF$SampleID<-gsub("SORT.","",SampleDF$Sample)
SampleDF$SampleID<-gsub(".bam","",SampleDF$SampleID)
pca_data2 <- merge(pca_data2, SampleDF, by.x = "Sample",by.y="SampleID")
pca_data2 <- pca_data2 %>% mutate(genetic_sex = case_when(avgBaf > 0.9 ~ "Female", TRUE ~ "Male"))
head(pca_data2)

## Plot the PCA after filtering
(PCAp2 <- ggplot(pca_data2, aes(x = PC1, y = PC2, colour = totalDP, shape = genetic_sex)) +
    theme_bw() + geom_point(alpha = 0.7, size = 3) +
    scale_colour_continuous(type = "viridis") )

ggsave("..//..//figures//sexmarkersPCA_filtered.pdf", plot = PCAp2, width = 5, height = 3.5)


#We can remove the 3 individuals that are in the middle of the PCA plot and then we are done with filtering and can analyse the data!
IDbadPCA<-pca_data2[pca_data2$PC1<1.2 & pca_data2$PC1>-1.2,]$avgBaf


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Estimate sex ratio----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#Get filtered dataset with sex and type and Area info
str(pca_data2)
str(popmap)
Sexed_inds<-merge(pca_data2,popmap, by.x="Sample", by.y="ID")#there is no need to filter since pca_data was already filtered for bad samples
Sexed_inds$Area<-as.factor(Sexed_inds$Area)
levels(Sexed_inds$Area)
#Sexed_inds$Area<- factor(Sexed_inds$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Ringstad'))
Sexed_inds$Area<- factor(Sexed_inds$Area, levels=c('Ringstad','Hitra','Austevoll', 'Arendal','Kristineberg'  ))
#Sexed_inds$type<-as.factor(Sexed_inds$type)
levels(Sexed_inds$type)
Sexed_inds$type2<-Sexed_inds$type
Sexed_inds$type2<-as.factor(Sexed_inds$type2)
levels(Sexed_inds$type2)<-c("eggs (T2)","female_control", "larvae (T3)" ,"male_control" ,"unsexed" )
#Sexed_inds$type2<-factor(Sexed_inds$type, levels=c("eggs (T2)","female_control", "larvae (T3)" ,"male_control" ,"unsexed" ))
  summary(Sexed_inds)

#define sex colors for plotting
sexcolors<-c("goldenrod1","mediumpurple1")

#Sample size
table(Sexed_inds$type)

## Plots----
#Overall sex ratio
ggplot(Sexed_inds,aes(x=genetic_sex,fill=genetic_sex))+
  geom_bar()+
  scale_fill_manual("Sex",values=sexcolors)+
  labs(title="Overall sex ratio", x="Genetic sex", y="Sample count")

  
#sex ratio by type and location, 
ggplot(Sexed_inds[Sexed_inds$type!="female_control" & Sexed_inds$type!="male_control" ,],aes(x=genetic_sex,fill=genetic_sex))+
  geom_bar()+
  scale_fill_manual("Sex",values=sexcolors)+
  #facet_wrap(fct_rev(Area)~type)
  facet_grid(cols=vars(type), rows=vars(Area))+
  labs(title="Sex ratio by type of sample and location", x="Genetic sex", y="Sample count")

#
#eggs and larvae alone first 
ggplot(Sexed_inds[Sexed_inds$type2=="eggs (T2)" | Sexed_inds$type2=="larvae (T3)"  ,],aes(x=genetic_sex,fill=genetic_sex))+
  geom_bar()+
  scale_fill_manual("Sex",values=sexcolors)+
  #facet_wrap(fct_rev(Area)~type)
  facet_grid(cols=vars(type2), rows=vars(Area))+
  labs(title="", x="Genetic sex", y="Sample count")
#eggs and larvae by nest
ggplot(Sexed_inds[Sexed_inds$type2=="eggs (T2)" | Sexed_inds$type2=="larvae (T3)"  ,],aes(x=genetic_sex,fill=genetic_sex))+
  geom_bar()+
  scale_fill_manual("Sex",values=sexcolors)+
  #facet_wrap(fct_rev(Area)~type)
  facet_grid(cols=vars(nest), rows=vars(Area))+
  labs(title="", x="Genetic sex", y="Sample count")


#unsexed by timepoint
ggplot(Sexed_inds[ Sexed_inds$type=="unsexed" ,],aes(x=genetic_sex,fill=genetic_sex))+
  geom_bar()+
  scale_fill_manual("Sex",values=sexcolors)+
  #facet_wrap(fct_rev(Area)~type)
  facet_grid(cols=vars(Timepoint), rows=vars(Area))+
  labs(title=" ", x="Genetic sex", y="Sample count")


#size of male and female unsexed individuals
ggplot(Sexed_inds[Sexed_inds$type=="unsexed",],aes(x=genetic_sex,y=size, fill=genetic_sex))+
  geom_boxplot()+
  geom_point(shape=21,alpha=0.4, position = position_jitterdodge())+
  scale_fill_manual("Sex",values=sexcolors)+
  labs(title="Body size of immature individuals", x="Genetic sex", y="Body size (mm)")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##Statistical test for the eggs/ larvae and unsexed sex ratios -----------------
str(Sexed_inds)
#I need to format the data for a binomial model!
#let's start with eggs/larvae in the same model
Sex_EL<-Sexed_inds[Sexed_inds$type=="larvae" | Sexed_inds$type=="eggs",]
T_EL<-as.data.frame(with(Sex_EL, table(Area,type,genetic_sex)))
F_EL<-T_EL[T_EL$genetic_sex=="Female",]
F_EL<-F_EL[,-3]
colnames(F_EL)[3]<-"Females"
M_EL<-T_EL[T_EL$genetic_sex=="Male",]
M_EL<-M_EL[,-3]
colnames(M_EL)[3]<-"Males"
B_EL<-merge(F_EL,M_EL,by=c("Area","type"))
B_EL<-B_EL[B_EL$Area!="Arendal" & B_EL$Area!="Austevoll",]
B_EL
#This new data formats allows a new type of plot that I find interesting
B_EL$N<-B_EL$Females+B_EL$Males
B_EL$pm<-B_EL$Males/B_EL$N
B_EL$CI<-1.96*sqrt(B_EL$pm*(1-B_EL$pm)/B_EL$N) #confidence interval using the normal approximation, justified by large sample size and high p
str(B_EL)
BEL<-B_EL[B_EL$type=="larvae" | B_EL$type=="eggs",]
levels(B_EL$type)<-c("Eggs (T2)","Larvae (T3)")
B_EL
sum(B_EL$N)

#make one at the nest level too
N_EL<-as.data.frame(with(Sex_EL, table(Area,nest,type,genetic_sex)))
N_EL<-N_EL[which(N_EL$Freq!=0),]
head(N_EL)
N_EL
F_EL<-N_EL[which(N_EL$genetic_sex=="Female"),]
F_EL<-F_EL[,-4]
colnames(F_EL)[4]<-"Females"
M_EL<-N_EL[which(N_EL$genetic_sex=="Male"),]
M_EL<-M_EL[,-4]
colnames(M_EL)[4]<-"Males"
Nest_EL<-merge(F_EL,M_EL,by=c("Area","nest","type"))
Nest_EL<-Nest_EL[which(Nest_EL$Area!="Arendal" & Nest_EL$Area!="Austevoll"),]
Nest_EL
#This new data formats allows a new type of plot that I find interesting
Nest_EL$N<-Nest_EL$Females+Nest_EL$Males
Nest_EL$pm<-Nest_EL$Males/Nest_EL$N
Nest_EL$CI<-1.96*sqrt(Nest_EL$pm*(1-Nest_EL$pm)/Nest_EL$N) #confidence interval using the normal approximation, justified by large sample size and high p
str(Nest_EL)
Nest_EL<-Nest_EL[which(Nest_EL$type=="larvae" | Nest_EL$type=="eggs"),]
levels(Nest_EL$type)<-c("Eggs (T2)","Larvae (T3)")
Nest_EL
sum(Nest_EL$N)
Nest_EL$Nest.unique<-paste0(Nest_EL$Area,"_N_",Nest_EL$nest)

# FIGURE 6 a----
## Plot egg/larvae SR----
A<-ggplot(B_EL, aes(y=pm, x=type, group=Area:type))+
  geom_point( size=4, aes(shape=type,fill=type))+
  ylim(0,1)+
  geom_hline(yintercept=0.5, alpha=0.5)+
  scale_shape_manual("Sample type",values=c(21,24))+
  scale_fill_manual("Sample type",values=c("royalblue",'orangered3'))+
  geom_errorbar(width=0.2,aes(ymin=pm-CI, ymax=pm+CI))+
  facet_wrap(~Area, ncol=3)+
  #theme(text=element_text(size=24))+
  theme_bw(base_size=13)+
  theme(legend.position = "bottom",plot.margin = unit(c(0.2,0.5,0.2,0.5), "cm"))+
  labs(y="Proportion of males", x="Sample type", title="")
#ggsave("..//..//..//Manuscripts//sex_det//figures//EggRatio.pdf", width = 8, height = 3.5)
A

## Plot egg/larvae SR by nest----
str(Nest_EL)
A2<-ggplot(Nest_EL, aes(y=pm, x=factor(nest), group=type))+
  geom_point( size=4, aes(shape=type,fill=type), position = position_dodge(width=0.2))+
  ylim(0,1)+
  geom_hline(yintercept=0.5, alpha=0.5)+
  scale_shape_manual("Sample type",values=c(21,24))+
  scale_fill_manual("Sample type",values=c("royalblue",'orangered3'))+
  geom_errorbar(width=0.2,aes(ymin=pm-CI, ymax=pm+CI),position = position_dodge(width=0.2))+
  facet_wrap(~Area, ncol=3,scales = "free_x")+
  #theme(text=element_text(size=24))+
  theme_bw(base_size=13)+
  theme(legend.position = "bottom",plot.margin = unit(c(0.2,0.5,0.2,0.5), "cm"))+
  labs(y="Proportion of males", x="Sample type", title="")
A2
ggsave("..//..//..//..//Manuscripts//sex_det//figures//EggRatio_nest.pdf", width = 8, height = 3.5)

## estimates and SE----
#Overall sex ratio?
modall<-glm(data=B_EL,family="binomial", cbind(Males,Females)~1)
summary(modall)
exp(-0.11980)/(1+exp(-0.11980)) #overall sex ratio 0.47
exp(-0.11980+0.08284  )/(1+exp(-0.11980+0.08284  ))  - exp(-0.11980)/(1+exp(-0.11980)) #SE
exp(-0.11980+0.08284 *qnorm(0.975) )/(1+exp(-0.11980+0.08284*qnorm(0.975)  )) #CI high
exp(-0.11980-0.08284 *qnorm(0.975) )/(1+exp(-0.11980-0.08284*qnorm(0.975)  )) #CI low
as.data.frame(effect("allEfects",modall,se=list(compute=T,level=0.95)))

mod1<-glm(data=B_EL,family="binomial", cbind(Males,Females)~Area*type)
Anova(mod1, type=3)
#non-sign interaction so...
mod2<-glm(data=B_EL,family="binomial", cbind(Males,Females)~Area+type)
summary(mod2)
exp(0.008517)/(1+exp(0.008517)) # intercept is 0.502 once translatd into p (Ringstad)
exp(0.008517+0.168578)/(1+exp(0.008517+0.168578)) #upper SE
exp(0.008517-0.168578)/(1+exp(0.008517-0.168578)) #lower
exp(0.008517)/(1+exp(0.008517)) - exp(0.008517-0.168578)/(1+exp(0.008517-0.168578))
exp(0.008517+0.168578)/(1+exp(0.008517+0.168578)) - exp(0.008517)/(1+exp(0.008517)) # I find SE =0.042
exp(0.008517-0.372553)/(1+exp(0.008517-0.372553)) #Kberg... 0.41 is it female biased??
Anova(mod2)

plot(effect("Area",mod2)) 
plot(effect("type",mod2))
as.data.frame(effect("Area",mod2,se=list(compute=T,level=0.95)))


coord<-read.delim("Coord.txt", header=T)
BEL_C<-merge(B_EL,coord,by="Area")

#no-signif interaction so removed
mod3<-glm(data=BEL_C,family="binomial", cbind(Males,Females)~Lat+type)
summary(mod3)
exp(-2.47753)/(1+exp(-2.47753))
Anova(mod3, type=2)
plot(effect("Lat",mod3))

# Unsexed adults ----
#Do same with unsexed adults
Sex_U<-Sexed_inds[Sexed_inds$type=="unsexed",]
T_U<-as.data.frame(with(Sex_U, table(Timepoint,genetic_sex)))
T_U$Timepoint<-as.character(T_U$Timepoint)
str(T_U)
F_U<-T_U[T_U$genetic_sex=="Female",]
F_U<-F_U[,-2]
colnames(F_U)[2]<-"Females"
M_U<-T_U[T_U$genetic_sex=="Male",]
M_U<-M_U[,-2]
colnames(M_U)[2]<-"Males"
B_U<-merge(F_U,M_U,by=c("Timepoint"))
B_U
#This new data formats allows a new type of plot that I find interesting
B_U$N<-B_U$Females+B_U$Males
B_U[4,]<-c(4,sum(B_U$Females),sum(B_U$Males),sum(B_U$N))
B_U[4,1]<-"ALL"
B_U$pm<-B_U$Males/B_U$N
B_U$CI<-1.96*sqrt(B_U$pm*(1-B_U$pm)/B_U$N) #confidence interval using the normal approximation, justified by large sample size and high p
str(B_U)
B_U

# FIGURE 6 b----
#plot by timepoint
B<-ggplot(B_U, aes(y=pm, x=Timepoint))+
  geom_point(data=B_U[B_U$Timepoint!="ALL",],aes( size=N),alpha=0.4, shape=22, fill="tan1")+
  geom_point(data=B_U[B_U$Timepoint=="ALL",],aes( size=N), shape=22, fill="tan1")+
  ylim(0,1)+
  scale_size_continuous("Sample size", range = c(2,5), breaks=c(30,96))+
  geom_hline(yintercept=0.5, alpha=0.5)+
 # scale_shape_manual("Sample type",values=c(21,24))+
  #scale_fill_manual("Sample type",values=c("royalblue",'orange'))+
  geom_errorbar(data=B_U[B_U$Timepoint!="ALL",],width=0.2,alpha=0.3,aes(ymin=pm-CI,ymax=pm+CI))+
  geom_errorbar(data=B_U[B_U$Timepoint=="ALL",],width=0.2,alpha=0.7,aes(ymin=pm-CI,ymax=pm+CI))+
  #facet_wrap(~Area, ncol=3)+
  theme_bw(base_size=13)+
  theme(legend.position = "bottom",plot.margin = unit(c(0.2,0.5,0.2,0.5), "cm"))+
  labs(y="Proportion of males", x="Time period", title="")
B
## merge plots----
#ggarrange(A,B, widths=c(2,1.), legend = "bottom")
cowplot::plot_grid(A,B,align="h", axis="bt", rel_widths = c(1.7,1), labels=c("a","b"))
ggsave("..//..//figures//SR_egg_unsexed.pdf", width = 9, height = 3.5)


#model
mod6<-glm(data=B_U,family="binomial", cbind(Males,Females)~Timepoint)
Anova(mod6, type=3)

str(SampleDFpretty)


#Appendix----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Code used to produce the sex specific genotypes from the origin ddRAD dataset
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# First I want to get the expected male and female genotypes from the ddRAD data!
sex_variants<-read.table("sex_candidates_final.txt", header=T) 
sex_variants
vcfRAD<-read.vcfR("maxmiss75_DP8_30.recode.vcf" )
RADvariants<-as.data.frame(vcfRAD@fix[, c("CHROM", "POS", "ID","REF", "ALT")])
head(RADvariants)
#Attach the genotype info
RADgenotypes <- extract.gt(vcfRAD)
RADvariants_g <- cbind(RADvariants, RADgenotypes)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#back to the original ddRAD data to keep only the variants that are also in the new data
RADvariants_sex_g<-RADvariants_g[ RADvariants_g$CHROM %in% variants$CHROM &
                                    RADvariants_g$POS %in% variants$POS , ]
RADvariants_sex_g
#based on this, I want to output the male and female genotypes for each locus.
#maybe the vcfR package has some helpful functions?
RADvariants_sex_g[,"SORT-01"] #this individual has only the ref allele, homozygous

#I need the sex info from a population map
maleID<-read.table("maleIDS.txt", header=F)
#I need to flip the RADvariant data vertically

RADvariants_sex_g2<-data.frame(colnames(RADvariants_sex_g)[6:183])
colnames(RADvariants_sex_g2)<-"ID"
#add sex info
RADvariants_sex_g2$sex<-"f"
RADvariants_sex_g2[RADvariants_sex_g2$ID %in% maleID$V1,]$sex<-"m"
RADvariants_sex_g2$sex<-as.factor(RADvariants_sex_g2$sex)
#add genotype of 9 loci
head(t(RADvariants_sex_g[1:9,6:183]))

rownames(RADvariants_sex_g2)<-colnames(RADvariants_sex_g)[6:183]

RADvariants_sex_g2<-cbind(RADvariants_sex_g2,t(RADvariants_sex_g[1:9,6:183]))
head(RADvariants_sex_g2)

RAD_sex_genotypes<-data.frame("CHR"=RADvariants_sex_g[,1],"POS"=RADvariants_sex_g[,2])
fem_Hom<-c()
mal_Het<-c()
fem_G<-c()
mal_G<-c()
for (i in 1:9) { 
  #count homozygous females and heterozygous males
  tab<-table(factor(RADvariants_sex_g2[,i+2]),RADvariants_sex_g2$sex)
  fHom<-tab[1] #0/0 count in f REF/REF
  fHet<-tab[2] #0/1 count in f
  mHom<-tab[3] #0/0 count in m
  mHet<-tab[4] #0/1 count in m REF/ALT
  fHomp<-fHom/(fHom+fHet)
  mHetp<-mHet/(mHet+mHom)
  fem_Hom<-c(fem_Hom,fHomp)
  mal_Het<-c(mal_Het,mHetp)
  # extract the reference genotype
  REF<-RADvariants_sex_g[i,4]
  ALT<-RADvariants_sex_g[i,5]
  fGeno<-paste(RADvariants_sex_g[i,4],RADvariants_sex_g[i,4], sep = "/")
  mGeno<-paste(RADvariants_sex_g[i,4],RADvariants_sex_g[i,5], sep = "/")
  fem_G<-c(fem_G,fGeno)
  mal_G<-c(mal_G,mGeno)
}

RAD_sex_genotypes$fGen<-fem_G
RAD_sex_genotypes$mGen<-mal_G
RAD_sex_genotypes$fem_Ho<-fem_Hom
RAD_sex_genotypes$mal_He<-mal_Het
RAD_sex_genotypes #now I know the main male and female genotypes for each of the 9 loci in the original ddRAD
#I can use this to classify individuals in the new dataset as males or females




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Code used during the investigation on the sample ID issue (-> discovered swapping of plates 7 & 8)



mapsize<-merge(maps2,sizes,by="ID")
mapsize$type<-"unknown"
realNTC<-as.character(c(36,72,108,215,251,322,394,501,537))
badNTC<-as.character(c(608,644,680,767,768))
susNTC<-as.character(c(179,358,465))
mapsize[mapsize$ID %in% realNTC,]$type<-"trusted NTC"
mapsize[mapsize$ID %in% badNTC,]$type<-"bad NTC (p7-p8)"
mapsize[mapsize$ID %in% susNTC,]$type<-"suspicious NTC"
mapsize$type<-as.factor(mapsize$type)
str(mapsize)
summary(mapsize)

ggplot(mapsize,aes(x=filesize, y=Percent_Mapped, fill=type))+
  geom_point(shape=21, size=2, alpha=0.6)+geom_text(data=mapsize[mapsize$type!="trusted NTC" &
                                                                   mapsize$filesize<3000 &
                                                                   mapsize$Percent_Mapped<75,],
                                                    aes(label=ID),hjust=-0.2, vjust=-0.3)



ggplot(data=mapsize[mapsize$Percent_Mapped>90  & mapsize$filesize<2000,],aes(x=filesize, y=Percent_Mapped, fill=type))+
  geom_point(shape=21, size=2, alpha=0.6)+geom_text(aes(label=ID),hjust=-0.2, vjust=-0.3)

#add to mapsize the within plate positional information

mapsize$plate<- ((as.numeric(mapsize$ID)-1)  %/% 96 )+1 #plate number
mapsize$plateCOL<- ((as.numeric(mapsize$ID)-1-96*(mapsize$plate-1))  %/% 8 ) +1 #column within plate
mapsize$plateROW<-  as.numeric(mapsize$ID)-96*(mapsize$plate-1)-8*(mapsize$plateCOL-1) #row within plate
mapsize$fsize.cut <- cut(mapsize$filesize, breaks = c(0,750,1500,3000,Inf), right = FALSE)
levels(mapsize$fsize.cut )<-c("0-750Kb","750Kb-1.5Mb","1.5Mb-3Mb",">3Mb")
mapsize$permap.cut <- cut(mapsize$Percent_Mapped, breaks = c(0,25,65,90,Inf), right = FALSE)
levels(mapsize$permap.cut )<-c("0-25%","25-65%","65-90%",">90%")

str(mapsize)
write_tsv(mapsize,"sampleinfo.txt")

## MAgic
ggplot(mapsize,aes(x=plateCOL,y=plateROW, fill=fsize.cut))+
  geom_point(shape=21,size=2,alpha=0.6,position = position_jitter(width=0.2, height=0.2))+
  scale_fill_manual("Fasta file size",values=c("brown2", "lightsalmon","#628cf3","olivedrab2"))+
  scale_y_reverse()+scale_x_continuous(breaks = seq(1,12,1))

ggplot(mapsize,aes(x=plateCOL,y=plateROW, fill=permap.cut))+
  geom_point(shape=21,size=2,alpha=0.6,position = position_jitter(width=0.2, height=0.2))+
  scale_fill_manual("% read mapped",values=c("brown2", "lightsalmon","#628cf3","olivedrab2"))+
  scale_y_reverse()+scale_x_continuous(breaks = seq(1,12,1))

ggplot(mapsize,aes(x=plateROW,y=filesize, group = plateROW))+
  geom_boxplot()+
  #scale_fill_manual("% read mapped",values=c("brown2", "lightsalmon","#628cf3","olivedrab2"))+
  scale_x_continuous(breaks = seq(1,12,1))

ggplot(mapsize,aes(x=plateCOL,y=filesize, group = plateCOL))+
  geom_boxplot()+
  #scale_fill_manual("% read mapped",values=c("brown2", "lightsalmon","#628cf3","olivedrab2"))+
  scale_x_continuous(breaks = seq(1,12,1))



ggplot(mapsize,aes(x=plateROW,y=Percent_Mapped, group = plateROW))+
  geom_boxplot()+
  #scale_fill_manual("% read mapped",values=c("brown2", "lightsalmon","#628cf3","olivedrab2"))+
  scale_x_continuous(breaks = seq(1,12,1))

ggplot(mapsize,aes(x=plateCOL,y=Percent_Mapped, group = plateCOL))+
  geom_boxplot()+
  #scale_fill_manual("% read mapped",values=c("brown2", "lightsalmon","#628cf3","olivedrab2"))+
  scale_x_continuous(breaks = seq(1,12,1))


