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
  'CurveFitFail|MultipleMeltPeak|NoAmplification', Raw_Data$Flags),]

# check Ct column has no NAs
Base_Data <- subset(Flags_Removed, Ct != "NA")

# create replicate and ID columns
Base_Data$Replicate = NA
Base_Data$ID = NA

# add replicate IDs
Base_Data <- mutate(Base_Data,
                   Replicate = case_when(
                     str_detect(Sample, "rep1") ~ "Rep1",
                     str_detect(Sample, "rep2") ~ "Rep2",
                     str_detect(Sample, "rep3") ~ "Rep3" ))

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

# pull "Assay" column out as a list.
Assay = dplyr::pull(Mean, Assay)

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


# test glm function on aminoglycoside?
Aminoglycoside_GLM <- glm(data = Aminoglycoside,
                          Fold_Change ~ Treatment_Stage)
summary(Aminoglycoside_GLM)

# number stages?
Gene_Expression <- Tidy_Data_All
Gene_Expression$Stage = NA
Gene_Expression <- mutate(Gene_Expression, Stage = case_when(
  str_detect(Treatment_Stage, "B_FeedPump") ~ "1",
  str_detect(Treatment_Stage, "C_Digester2") ~ "2",
  str_detect(Treatment_Stage, "D_Digester4") ~ "3" ,
  str_detect(Treatment_Stage, "E_CentrifugeIn") ~ "4",
  str_detect(Treatment_Stage, "F_NewDrySolid") ~ "5",
  str_detect(Treatment_Stage, "G_TWDrySolid") ~ "6",
  str_detect(Treatment_Stage, "H_FSTDist") ~ "7"))

