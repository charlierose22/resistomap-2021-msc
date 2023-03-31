# import raw dataset
rawdata <- readr::read_csv("Data-Raw/rawdata.csv") %>% janitor::clean_names()

# Remove unsatisfactory flags.
flags_removed <- rawdata[!grepl('CurveFitFail|MultipleMeltPeak|NoAmplification',
                                 rawdata$flags),]
# Remove flags column, as it's no longer needed.
flags_removed$flags = NULL

# Add individual sample locations for removing "solo" results.
# "Solo" results are when out of the three repeats, there was only amplification in one sample.
flags_removed$id = NA
flags_removed$replicate = NA

# Fill in replicate IDs.
flags_removed <- mutate(flags_removed,
                        replicate = case_when(
                          str_detect(sample, "rep1") ~ "1",
                          str_detect(sample, "rep2") ~ "2",
                          str_detect(sample, "rep3") ~ "3" ))

# Fill in sample IDs.
# I swapped ASP Distribution to be Sample 1 as it's the earliest stage.
flags_removed <- mutate(flags_removed, id = case_when(
  str_detect(sample, "1A") ~ "A",
  str_detect(sample, "2A") ~ "B",
  str_detect(sample, "3A") ~ "C" ,
  str_detect(sample, "4A") ~ "D",
  str_detect(sample, "5A") ~ "E",
  str_detect(sample, "6A") ~ "F",
  str_detect(sample, "7A") ~ "ASP",
  str_detect(sample, "8A") ~ "G"))

# Remove "solo" results.
solo_removed <- ddply(flags_removed, c("assay", "id"),
                      function(d) {if (nrow(d) > 1) d else NULL})

# Remove Tm and Efficiency columns as we're focusing on cycle threshold.
solo_removed$tm = NULL
solo_removed$efficiency = NULL

# We now have to remove entries with no amplification at the ASP location.
# ASP is the first location during the treatment process we have data for.
# Make a list of unique genes that have data for the ASP location, then filter by that. list.
list_assay_id <- unique(solo_removed[grepl("ASP",solo_removed$id),]$assay)
only_asp <- dplyr::filter(solo_removed, solo_removed$assay %in% list_assay_id)

# Mutate the sample column to not begin with a number, for easier coding and recognition.
only_asp$sample_id = NA
mutated_sample_id <- mutate(only_asp, sample_id = case_when(
  str_detect(sample, "1A-rep1") ~ "A_rep1",
  str_detect(sample, "1A-rep2") ~ "A_rep2",
  str_detect(sample, "1A-rep3") ~ "A_rep3",
  str_detect(sample, "2A-rep1") ~ "B_rep1",
  str_detect(sample, "2A-rep2") ~ "B_rep2",
  str_detect(sample, "2A-rep3") ~ "B_rep3",
  str_detect(sample, "3A-rep1") ~ "C_rep1" ,
  str_detect(sample, "3A-rep2") ~ "C_rep2" ,
  str_detect(sample, "3A-rep3") ~ "C_rep3" ,
  str_detect(sample, "4A-rep1") ~ "D_rep1",
  str_detect(sample, "4A-rep2") ~ "D_rep2",
  str_detect(sample, "4A-rep3") ~ "D_rep3",
  str_detect(sample, "5A-rep1") ~ "E_rep1",
  str_detect(sample, "5A-rep2") ~ "E_rep2",
  str_detect(sample, "5A-rep3") ~ "E_rep3",
  str_detect(sample, "6A-rep1") ~ "F_rep1",
  str_detect(sample, "6A-rep2") ~ "F_rep2",
  str_detect(sample, "6A-rep3") ~ "F_rep3",
  str_detect(sample, "7A-rep1") ~ "ASP_rep1",
  str_detect(sample, "7A-rep2") ~ "ASP_rep2",
  str_detect(sample, "7A-rep3") ~ "ASP_rep3",
  str_detect(sample, "8A-rep1") ~ "G_rep1",
  str_detect(sample, "8A-rep2") ~ "G_rep2",
  str_detect(sample, "8A-rep3") ~ "G_rep3"))

# Remove the original Sample ID column, as it has been replaced.
mutated_sample_id$sample = NULL

# Pivot the table wider.
mutated_wide <- mutated_sample_id %>%
  pivot_wider(names_from = assay, values_from = ct)

# Remove extra unneeded columns
mutated_wide$id = NULL
mutated_wide$replicate = NULL

# Create sample ID list.
sample_id <- mutated_wide$sample_id

# Pivot the table back to longer.
mutated_long <- mutated_wide %>%
  pivot_longer(
    cols = AY1:AY96, names_to = "assay", values_to = "ct")

# Pivot the table again, so Sample ID is across the top, rather than Assay.
mutated_wide_2 <- mutated_long %>%
  pivot_wider(
    names_from = sample_id, values_from = ct)

# Transpose the table.
transposed_ct <- t(mutated_wide_2[2:25])

# Make sure the top row is the Assay codes, by first extracting as a list.
assay_names <- t(mutated_wide_2[1])
colnames(transposed_ct) <- as.character(assay_names[1,])

# Calculate Delta Ct and make sure the output is as a data frame.
as.data.frame(transposed_ct)
delta_ct <- transposed_ct[ , 2:137] - transposed_ct[ , "AY1"]
df_delta_ct <- as.data.frame(delta_ct)

# Start with ASP values on their own.
delta_ct_asp <- head(df_delta_ct, 3)

# Turn the rownames into the first column to preserve them.
delta_ct_asp_rownames <- rownames_to_column(delta_ct_asp, "sample_id")

# Calculate the sum of each column.
delta_ct_asp_sum <- as.data.frame(colSums(delta_ct_asp_rownames[ , -1], na.rm = TRUE))

# Rename the resulting sum column.
delta_ct_asp_sum_rownames <- rownames_to_column(delta_ct_asp_sum, "assay")
names(delta_ct_asp_sum_rownames)[2] <- "sum_delta_ct"

# Pivot wider ready for mean calculation.
delta_ct_asp_sum_wide <- pivot_wider(delta_ct_asp_sum_rownames,
                                     names_from = "assay",
                                     values_from = "sum_delta_ct")

# Pivot longer for counting data entries.
delta_ct_asp_long <- pivot_longer(delta_ct_asp_sum_wide,
                                  cols = AY10:AY96,
                                  names_to = "assay",
                                  values_to = "delta_ct",
                                  values_drop_na = TRUE)

# Count the frequency for each gene.
delta_ct_asp_count <- as.data.frame(table(delta_ct_asp_long$assay))

# Pivot wider again.
delta_ct_asp_count_wide <- as.data.frame(pivot_wider(delta_ct_asp_count,
                                                     names_from = "Var1",
                                                     values_from = "Freq"))

# Calculate ASP means for Delta Delta Ct calculation.
delta_ct_asp_mean <- as.data.frame(delta_ct_asp_sum_wide * delta_ct_asp_count_wide)

# Calculate Delta Delta Ct.
delta_delta_ct <- df_delta_ct - as.list(delta_ct_asp_mean)

# Remove ASP values from main table.
delta_delta_ct_no_asp <- tail(delta_delta_ct, -3)

# Calculate 2^DDC, or is it better to log2 later?
ddct_power <- 2^-(delta_delta_ct_no_asp)

# Convert row names to a column of their own to protect them.
ddct_power_location <- rownames_to_column(ddct_power, "sample_id")

# Add replicate and ID columns.
ddct_p <- ddct_power_location
ddct_p$replicate = NA
ddct_p$treatment_stage = NA

# Rearrange the columns to make it easier.
rearranged_ddct_p <- subset(ddct_p, select = c(
  sample_id, treatment_stage, replicate, AY10:AY96))

rearranged_ddct_p <- mutate(rearranged_ddct_p,
                            replicate = case_when(
                              str_detect(sample_id, "rep1") ~ "rep1",
                              str_detect(sample_id, "rep2") ~ "rep2",
                              str_detect(sample_id, "rep3") ~ "rep3"))

# Add sample IDs in their column.
rearranged_ddct_p <- mutate(rearranged_ddct_p, treatment_stage = case_when(
  str_detect(sample_id, "A") ~ "A",
  str_detect(sample_id, "B") ~ "B",
  str_detect(sample_id, "C") ~ "C",
  str_detect(sample_id, "D") ~ "D",
  str_detect(sample_id, "E") ~ "E",
  str_detect(sample_id, "F") ~ "F",
  str_detect(sample_id, "G") ~ "G"))

# Remove the joined sample ID column, as it has been split.
rearranged_ddct_p$sample_id = NULL

# Pivot longer.
ddct_p_long <- pivot_longer(rearranged_ddct_p,
                            cols = AY10:AY96,
                            names_to = "assay",
                            values_to = "ddct_two_p",
                            values_drop_na = TRUE)

ddct_asp_power <- 2^-(delta_ct_asp_mean)
ddct_asp_mean_long <- pivot_longer(ddct_asp_power,
                                   cols = AY10:AY96,
                                   names_to = "assay",
                                   values_to = "asp_dct")

# Join the ASP and long table together.
ddct_p_with_asp <- full_join(ddct_p_long, ddct_asp_mean_long,
                            by = c("assay" = "assay"))



# Import assay information.
assay_information <- readr::read_csv("Data-Raw/assayinformation.csv") %>% 
  janitor::clean_names()

# Remove the columns not needed.
assay_information$forward_primer = NULL
assay_information$reverse_primer = NULL

# Change column names for easier code.
colnames(assay_information) = c("assay", "gene", "target_antibiotic")

# Remove Taxonomic genes for now, focus on resistance genes.
assay_information <- assay_information[!grepl('Taxanomic',
                                              assay_information$target_antibiotic),]

# Join the main table with the assay information.
annotated_ddct_p <- full_join(ddct_p_with_asp, assay_information,
                             by = c("assay" = "assay"))

# Remove any NAs
Annotated_DDCtP <- Annotated_DDCtP %>% drop_na()

# Rearrange columns.
Annotated_DDCtP <- subset(Annotated_DDCtP, select = c(
  Assay, Treatment_Stage, Gene, DDCtTwoP, ASP_DCt, Target_Antibiotic, Replicate))

# Calculate the standard errors.
DDCtP_Summary <- Annotated_DDCtP %>%
  group_by(Gene, Treatment_Stage) %>%
  summarise(mean = mean(DDCtTwoP),
            std = sd(DDCtTwoP),
            n = length(DDCtTwoP),
            se = std/sqrt(n))

Summary_DDCtP <- Annotated_DDCtP
Summary_DDCtP$Assay = NULL

# Potentially create a wide format table with the summary information.
DDCtP_Wide <- pivot_wider(Annotated_DDCtP,
                          names_from = Gene,
                          values_from = DDCtTwoP)

# Produce initial plot.
ggplot(data = Annotated_DDCtP, mapping =
         aes(Treatment_Stage, DDCtTwoP)) +
  geom_violin() +
  facet_wrap_paginate(facets = vars(Target_Antibiotic),
                      ncol = 4,
                      scales = "free_x")

# Split based on Treatment Stage.
Split <- split(Annotated_DDCtP, Annotated_DDCtP$Treatment_Stage)
A <- Split$A
B <- Split$B
C <- Split$C
D <- Split$D
E <- Split$E
F <- Split$F
G <- Split$G

# Rename dataset for model.
DDCtP_LM <- DDCtP_Wide
DDCtP_LM <- na.omit(DDCtP_LM)
str(DDCtP_LM)

# Create model.
# How does DDCtTwoP change across each location for each gene?
A_lm <- lm(DDCtTwoP ~ Gene, data = A)

# Use tidy to neaten table and add confidence intervals.
A_out_conf <- tidy(A_lm, conf.int = TRUE)

# Needed to strip term names.
# Strip out prefix in term column.
A_out_conf$nicelabs <- prefix_strip(A_out_conf$term, "Gene")

# Augment generates Cook's distance, residual values and fitted values.
A_out_aug <- augment(A_lm)

# Create new column to show location.
A_out_aug$Treatment_Stage = "A"

# Summary model values.
glance(A_lm)

# repeat for B.
B_lm <- lm(DDCtTwoP ~ Gene, data = B)
B_out_conf <- tidy(B_lm, conf.int = TRUE)
B_out_conf$nicelabs <- prefix_strip(B_out_conf$term, "Gene")
B_out_aug <- augment(B_lm)
B_out_aug$Treatment_Stage = "B"
glance(B_lm)

# repeat for C.
C_lm <- lm(DDCtTwoP ~ Gene, data = C)
C_out_conf <- tidy(C_lm, conf.int = TRUE)
C_out_conf$nicelabs <- prefix_strip(C_out_conf$term, "Gene")
C_out_aug <- augment(C_lm)
C_out_aug$Treatment_Stage = "C"
glance(C_lm)

# repeat for D.
D_lm <- lm(DDCtTwoP ~ Gene, data = D)
D_out_conf <- tidy(D_lm, conf.int = TRUE)
D_out_conf$nicelabs <- prefix_strip(D_out_conf$term, "Gene")
D_out_aug <- augment(D_lm)
D_out_aug$Treatment_Stage = "D"
glance(D_lm)

# repeat for E.
E_lm <- lm(DDCtTwoP ~ Gene, data = E)
E_out_conf <- tidy(E_lm, conf.int = TRUE)
E_out_conf$nicelabs <- prefix_strip(E_out_conf$term, "Gene")
E_out_aug <- augment(E_lm)
E_out_aug$Treatment_Stage = "E"
glance(E_lm)

# repeat for F.
F_lm <- lm(DDCtTwoP ~ Gene, data = F)
F_out_conf <- tidy(F_lm, conf.int = TRUE)
F_out_conf$nicelabs <- prefix_strip(F_out_conf$term, "Gene")
F_out_aug <- augment(F_lm)
F_out_aug$Treatment_Stage = "F"
glance(F_lm)

# repeat for G.
G_lm <- lm(DDCtTwoP ~ Gene, data = G)
G_out_conf <- tidy(G_lm, conf.int = TRUE)
G_out_conf$nicelabs <- prefix_strip(G_out_conf$term, "Gene")
G_out_aug <- augment(G_lm)
G_out_aug$Treatment_Stage = "G"
glance(G_lm)

# Join all tables together.
Full_lm_out <- rbind.data.frame(A_out_aug,
                                B_out_aug,
                                C_out_aug,
                                D_out_aug,
                                E_out_aug,
                                F_out_aug,
                                G_out_aug)

# Join with the assay information.
Annotated_Full_Model <- full_join(Full_lm_out, Assay_Information,
                                  by = c("Gene" = "Gene"))

Annotated_Full_Model <- Annotated_Full_Model %>% drop_na()
Annotated_Full_Model$Assay = NULL

# Re-do summary table.
Annotated_Full_Summary <- Annotated_Full_Model %>%
  group_by(Gene, Treatment_Stage) %>%
  summarise(mean = mean(DDCtTwoP),
            std = sd(DDCtTwoP),
            n = length(DDCtTwoP),
            se = std/sqrt(n))

# Join summary and model table.
All_Data <- full_join(Annotated_Full_Summary, Annotated_Full_Model,
                      by = c("Gene" = "Gene", "Treatment_Stage" = "Treatment_Stage"))

# create a table for tidy data
write.csv(All_Data, "AllData.csv", row.names = FALSE)
# ------------GRAPHS------------------
# Loop in ggplot.
n = 150
pdf("resistomap.pdf", paper= 'A4r', width = 8, height = 6)
for (i in 1:n) {
  print(ggplot(data = All_Data,
               aes(x = Treatment_Stage,
                   y = mean)) +
          geom_point() +
          facet_wrap_paginate(vars(Gene), scales = "free",
                              ncol = 1, nrow = 1, page = i) +
          geom_errorbar(aes(x = Treatment_Stage,
                            ymin = mean,
                            ymax = mean),
                        width = .3) +
          geom_errorbar(aes(x = Treatment_Stage,
                            ymin = mean - se,
                            ymax = mean + se),
                        width = .5) +
          labs(x = "Treatment Stage",
               y = "Normalised Gene Expression") +
          theme_bw(base_size = 12) + facet_wrap_paginate(vars(Gene), scales = "free",
                                                         ncol = 1, nrow = 1, page = i))}
dev.off()

# heatmap
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
ggplot(All_Data, aes(y = Gene, x = Treatment_Stage, fill = mean)) +
  geom_tile(aes(fill = mean)) +
  scale_fill_gradientn(colours = mycol, trans = "log") +
  labs(x = "Wastewater Treatment Stage", 
       y = "Gene", 
       fill = "Gene Prevalence") +
  theme_bw()
ggsave("Figure/HeatmapTotalGene.png", width = 6, height = 30)

#split table based on target antibiotic??
Split <- split(All_Data, All_Data$Target_Antibiotic)
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
Tetracycline <- Split$Tetracycline
Trimethoprim <- Split$Trimethoprim

# seperate heatmaps?
