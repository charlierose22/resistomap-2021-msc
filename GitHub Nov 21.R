#----------------SETUP-------------------

setwd("~/GitHub/DNA-Analysis-2021")
library(plyr)
library(tidyverse)
library(readr)
library(ggforce)
library(colorBlindness)

# import raw dataset
Raw_Data <- read_csv("Data-Raw/rawdata.csv")

# remove unsatisfactory flags
Flags_Removed <- Raw_Data[!grepl(
  'CurveFitFail|MultipleMeltPeak|NoAmplification', Raw_Data$Flags), ]

# ---------------------INDIVIDUALS---------------------------

# rename for stats without means.
Individual <- Flags_Removed
Indv_Solo_Removed$Flags = NULL

# Add individual sample locations for removing solo results
Individual$ID = NA
Individual$Replicate = NA

# add replicate IDs
Individual <- mutate(Individual,
                     Replicate = case_when(
                       str_detect(Sample, "rep1") ~ "1",
                       str_detect(Sample, "rep2") ~ "2",
                       str_detect(Sample, "rep3") ~ "3" ))

# add sample IDs (swapped ASP Distribution to be Sample 1 as it's the earliest stage)
Individual <- mutate(Individual, ID = case_when(
  str_detect(Sample, "1A") ~ "B_FeedPump",
  str_detect(Sample, "2A") ~ "C_Digester2",
  str_detect(Sample, "3A") ~ "D_Digester4" ,
  str_detect(Sample, "4A") ~ "E_CentrifugeIn",
  str_detect(Sample, "5A") ~ "F_NewDrySolid",
  str_detect(Sample, "6A") ~ "G_TWDrySolid",
  str_detect(Sample, "7A") ~ "A_ASPDist",
  str_detect(Sample, "8A") ~ "H_FSTDist"))

# remove entries with data for only one replicate
Indv_Solo_Removed <- ddply(Individual, c("Assay", "ID"),
                           function(d) {if (nrow(d) > 1) d else NULL})

# remove entries with no data for ASP
Indv_Wider_ASP <- pivot_wider(Indv_Solo_Removed, )



# remove Tm and Efficiency columns
Indv_Solo_Removed$Tm = NULL
Indv_Solo_Removed$Efficiency = NULL

# change sample column to have locations in and not begin with a number.
Indv_Solo_Removed$SampleID = NA
Individual_Mutated <- mutate(Indv_Solo_Removed, SampleID = case_when(
  str_detect(Sample, "1A-rep1") ~ "B1_FeedPump",
  str_detect(Sample, "1A-rep2") ~ "B2_FeedPump",
  str_detect(Sample, "1A-rep3") ~ "B3_FeedPump",
  str_detect(Sample, "2A-rep1") ~ "C1_Digester2",
  str_detect(Sample, "2A-rep2") ~ "C2_Digester2",
  str_detect(Sample, "2A-rep3") ~ "C3_Digester2",
  str_detect(Sample, "3A-rep1") ~ "D1_Digester4" ,
  str_detect(Sample, "3A-rep2") ~ "D2_Digester4" ,
  str_detect(Sample, "3A-rep3") ~ "D3_Digester4" ,
  str_detect(Sample, "4A-rep1") ~ "E1_CentrifugeIn",
  str_detect(Sample, "4A-rep2") ~ "E2_CentrifugeIn",
  str_detect(Sample, "4A-rep3") ~ "E3_CentrifugeIn",
  str_detect(Sample, "5A-rep1") ~ "F1_NewDrySolid",
  str_detect(Sample, "5A-rep2") ~ "F2_NewDrySolid",
  str_detect(Sample, "5A-rep3") ~ "F3_NewDrySolid",
  str_detect(Sample, "6A-rep1") ~ "G1_TWDrySolid",
  str_detect(Sample, "6A-rep2") ~ "G2_TWDrySolid",
  str_detect(Sample, "6A-rep3") ~ "G3_TWDrySolid",
  str_detect(Sample, "7A-rep1") ~ "A1_ASPDist",
  str_detect(Sample, "7A-rep2") ~ "A2_ASPDist",
  str_detect(Sample, "7A-rep3") ~ "A3_ASPDist",
  str_detect(Sample, "8A-rep1") ~ "H1_FSTDist",
  str_detect(Sample, "8A-rep2") ~ "H2_FSTDist",
  str_detect(Sample, "8A-rep3") ~ "H3_FSTDist"))

Individual_Mutated$Sample = NULL

# Pivot table
Indv_Wide <- Individual_Mutated %>%
  pivot_wider(names_from = Assay, values_from = Ct)

# Create sample ID table
SampleID <- Indv_Wide$SampleID

# Calculate DCt
#Indv_DCt <- Indv_Wide %>%
#mutate(across(c(AY10:AY96),
#   ~.x - "AY1",
#   .names = "{.col}-AY1"),)
# <- Indv_DCt[,grep("-AY1|Sample", colnames(Indv_DCt))]

# calculate DDCt

# pivot back to long
Indv_Long <- Indv_Wide %>% pivot_longer(
  cols = AY1:AY96, names_to = "Assay", values_to = "DeltaCt")

# import assay information
Assay_Information <- read_csv("Data-Raw/assayinformation.csv")

# remove columns not needed
Assay_Information$`Forward Primer` = NULL
Assay_Information$`Reverse Primer` = NULL

# change column names for easier code
colnames(Assay_Information) = c("Assay", "Gene", "Target_Antibiotic")

# join assay information
Annotated_Indv <- full_join(Indv_Long, Assay_Information,
                            by = c("Assay" = "Assay"))
Tidy_Annotated_Indv <- Annotated_Indv %>% drop_na()



# ---------------------MEANS---------------------------

# check Ct column has no NAs
Base_Data <- subset(Flags_Removed, Ct != "NA")
Base_Data$Flags = NULL

# create replicate and ID columns
Base_Data$Replicate = NA
Base_Data$ID = NA

# add replicate IDs
Base_Data <- mutate(Base_Data,
                    Replicate = case_when(
                      str_detect(Sample, "rep1") ~ "A",
                      str_detect(Sample, "rep2") ~ "B",
                      str_detect(Sample, "rep3") ~ "C" ))

# add sample IDs (swapped ASP Distribution to be Sample 1 as it's the earliest stage)
Base_Data <- mutate(Base_Data, ID = case_when(
  str_detect(Sample, "1A") ~ "B_FeedPump",
  str_detect(Sample, "2A") ~ "C_Digester2",
  str_detect(Sample, "3A") ~ "D_Digester4" ,
  str_detect(Sample, "4A") ~ "E_CentrifugeIn",
  str_detect(Sample, "5A") ~ "F_NewDrySolid",
  str_detect(Sample, "6A") ~ "G_TWDrySolid",
  str_detect(Sample, "7A") ~ "A_ASPDist",
  str_detect(Sample, "8A") ~ "H_FSTDist"))

# remove entries with data for only one replicate
Solo_Removed <- ddply(Base_Data, c("Assay", "ID"),
                      function(d) {if (nrow(d) > 1) d else NULL})

# calculate means
Mean <- Solo_Removed %>%
  group_by(Assay, ID) %>%
  summarise_at(vars(Ct, Tm, Efficiency),
               list(Mean = mean))

# remove Tm and Efficiency columns
Ct_Only <- Mean
Ct_Only$Tm_Mean = NULL
Ct_Only$Efficiency_Mean = NULL

# remove assays which don't have a value for ASP
Ct_Only_Wide <- Ct_Only %>%
  pivot_wider(names_from = ID, values_from = Ct_Mean)

Ct_Only_ASP = subset(Ct_Only_Wide, Ct_Only_Wide$A_ASPDist != "NA")

ASP_Only = Ct_Only_ASP %>% pivot_longer(
  cols = A_ASPDist:H_FSTDist, names_to = "ID", values_to = "Ct_mean")

# create a wide table and transpose
Wide_Mean_Ct <- ASP_Only %>%
  pivot_wider(names_from = ID, values_from = Ct_mean)

Transposed_Wide_Mean_Ct <- t(Wide_Mean_Ct[2:9])

Transposed_Wide_Mean_Ct[is.na(Transposed_Wide_Mean_Ct)] <- 0
# keeps numeric but we've lost assay codes

# make sure the top row is the Assay codes
Assay_Names <- t(Wide_Mean_Ct[1])
colnames(Transposed_Wide_Mean_Ct) <- as.character(Assay_Names[1,])

# change 0 to NA
Transposed_Wide_Mean_Ct[Transposed_Wide_Mean_Ct == 0] <- NA

# Calculate DCt and DDCt values.
Delta_Ct <- Transposed_Wide_Mean_Ct[ , 2:137] - Transposed_Wide_Mean_Ct[ , "AY1"]
DF_Delta_Ct <- as.data.frame(Delta_Ct)
Delta_Delta_Ct <- DF_Delta_Ct - as.list(DF_Delta_Ct[1, ])
Delta_Delta_Ct_No_ASP <- Delta_Delta_Ct[-1, ]
DDCt_Power <- 2^-(Delta_Delta_Ct_No_ASP)

# add ID variable back to row names.
ID_DDCt_Power <- rownames_to_column(DDCt_Power, var = "ID")
Assay_DDCt_Power <- column_to_rownames(ID_DDCt_Power, "ID")

# transpose table, to add assay annotations.
DDCt_Power_Long <- as.data.frame(t(Assay_DDCt_Power))

# import assay information
Assay_Information <- read_csv("Data-Raw/assayinformation.csv")

# remove columns not needed
Assay_Information$`Forward Primer` = NULL
Assay_Information$`Reverse Primer` = NULL

# change column names for easier code
colnames(Assay_Information) = c("Assay", "Gene", "Target_Antibiotic")

# join assay information with DDCt values.
Assays_DDCt_Power_Long <- rownames_to_column(DDCt_Power_Long, var = "Assay")
Annotated_DDCt <- full_join(Assays_DDCt_Power_Long, Assay_Information,
                            by = c("Assay" = "Assay"))
Tidy_Annotated_DDCt <- Annotated_DDCt %>% drop_na()

# pivot table longer to create a row of sample IDs
Tidy_Data_All <- Tidy_Annotated_DDCt %>%
  pivot_longer(cols = B_FeedPump:H_FSTDist,
               names_to = "Treatment_Stage",
               values_to = "Fold_Change")
group_by(Tidy_Data_All, Gene, Target_Antibiotic)

# counts of genes? All 7?
Genes_Found <- as.matrix(Tidy_Data_All$Gene)
Genes_Found <- t(Genes_Found)
Genes_Found <- Genes_Found[!duplicated(Genes_Found),]
Genes_Found <- as.data.frame(Genes_Found)
aggregate(Genes_Found, by = Genes_Found, length) [1:ncol(Genes_Found) + 1]

# ---------------------GRAPHS-----------------------------

# produce initial graph based on tidy data
pdf("Figures/All_Fold_Expressions.pdf", width = 10, height = 50)
ggplot(data = Tidy_Data_All, mapping = aes(x = Gene, y = Fold_Change)) +
  geom_point(aes(color = Treatment_Stage)) +
  facet_wrap_paginate(facets = vars(Target_Antibiotic),
                      ncol = 1,
                      scales = "free_x") +
  labs(x = "Gene",
       y = "Relative Fold Expression",
       color = "Treatment Stage") +
  theme_classic()
dev.off()

#split based on target antibiotic
Split <- split(Tidy_Data_All, Tidy_Data_All$Target_Antibiotic)
Aminoglycoside <- Split$Aminoglycoside
Beta_Lactam <- Split$`Beta Lactam`
Integrons <- Split$Integrons
MDR <- Split$MDR
MGE <- Split$MGE
MLSB <- Split$MLSB
Other <- Split$Other
Phenicol <- Split$Phenicol
Quinolone <- Split$Quinolone
Sulfonamide <- Split$Sulfonamide
Taxonomic <- Split$Taxanomic
Tetracycline <- Split$Tetracycline
Trimethoprim <- Split$Trimethoprim

# Focus on Aminoglycoside points first.
ggplot(data = Aminoglycoside,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/Aminoglycoside_1.png", width = 12, height = 8)

# Copy line graohs for each of the Targets
ggplot(data = Beta_Lactam,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/Beta_Lactam_1.png", width = 12, height = 8)

ggplot(data = Integrons,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/Integrons_1.png", width = 12, height = 8)

ggplot(data = MGE,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/MGE_1.png", width = 12, height = 8)

ggplot(data = MLSB,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/MLSB_1.png", width = 12, height = 8)

ggplot(data = Other,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/Other_1.png", width = 12, height = 8)

ggplot(data = Phenicol,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/Phenicol_1.png", width = 12, height = 8)

ggplot(data = Quinolone,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/Quinolone_1.png", width = 12, height = 8)

ggplot(data = Sulfonamide,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/Sulfonamide_1.png", width = 12, height = 8)

ggplot(data = Taxonomic,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/Taxonomic_1.png", width = 12, height = 8)

ggplot(data = Tetracycline,
       mapping = aes(x = Treatment_Stage, y = Fold_Change,
                     group = Gene, color = as.factor(Gene))) +
  geom_point(alpha = 0.5, aes(color = Gene)) +
  geom_line(aes(color =  Gene)) +
  labs(x = "Treatment Stage",
       y = "Relative Fold Expression",
       colour = "Gene") +
  theme_light()
dev.off()
ggsave("Figures/Tetracycline_1.png", width = 12, height = 8)
