# 1. Basic set up & loading data ----

# clean environment
rm(list=ls())

# load libraries
library(ggplot2)
library(dplyr)
library(multcompView)
library(agricolae)
library(ggsignif)
library(ggpubr)
library("report")
library("readxl")

# Set workingdirectory
setwd("/home/melanie/working_directory/D36E_assays/ATTA_spinach")

# generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")

# filter for samples of interest 
pst_raw <- read_xlsx("/home/melanie/working_directory/D36E_assays/CFU_data/LDM20230330_CFU_count_calculations_AGL1_effect_on_pst.xlsx")

# 2. Select/ filter data ----
# 2.1 Bacteria strains as factors
pst_raw$Bacteria.strain <- as.factor(pst_raw$`Bacterial strain`)
pst_raw$Combination <- as.factor(pst_raw$Combination)
#pst_raw$Date <- as.factor(pst_raw$dpi)
pst_raw$Date <- as.factor(
  ifelse(pst_raw$dpi == -2, "AminusTwo",
  ifelse(pst_raw$dpi == 0, "Bzero",
  ifelse(pst_raw$dpi == 3, "Cthree", pst_raw$dpi)))
)

# remove rows with NA (if extra rows are loaded)
# pst_raw <- pst_raw[1:210,]

# 2.2 log10 transformation of data for better presentation

#MIN <- min(pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`[pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`!=0]) # set lowest CFU counts as 0
pst_raw$CFU.log10 <- log10(pst_raw$`CFU/ml/ 1cm^2`+1)



# 2.3 make df containing only data I want to plot: Pst strain, CFU log10 transormed
Pst_strain <- as.factor(pst_raw$Bacteria.strain)
combination <- pst_raw$Combination
Pst_CFU <- pst_raw$CFU.log10
Pst_exp_date <- as.factor(pst_raw$Date)
#Pst_OD <- pst_raw$OD6000.infiltrated.0dpi
#Pst_comment <- pst_raw$Comment

simple_df <- data.frame(Pst_strain, combination, Pst_exp_date, Pst_CFU)

# 2.4 Filter for group of interest
simple_df_agl <- simple_df[simple_df$Pst_strain == "AGL1",]
simple_df_agl <- simple_df_agl[simple_df_agl$Pst_CFU != 0, ]

co <- c("#5b5d5cff", "#44AA99", "#E69F00")  # Define colors

p_agl <- ggplot(simple_df_agl, aes(x = as.factor(combination), y = Pst_CFU, fill = Pst_exp_date)) +
  geom_boxplot(outlier.colour = "red", alpha = 0.9, 
               color = "black",  # Black outlines, whiskers, and median
               fatten = 2) +  # Thicker median line
  labs(x = "Pst_exp_date", y = "log10 (CFU/cm2)") +
  geom_point(aes(color = Pst_exp_date), position = position_jitterdodge(0.1)) +  # Scatter colored by date
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_fill_manual(values = c("white", "white", "white")) +  # Fill colors for boxes
  scale_color_manual(values = co) +  # Colors for scatter points
  ylim(0, 8)

p_agl

ggsave(filename =  paste0("ATTA_CFU_count_Agl_boxplot.svg"), 
       plot = p_agl,
       device = "svg")


# 4. Statistical analysis across samples using ANOVA+TUKEY ----
# 4.1. ANOVA to analyse whether there are differences between the groups
#      H0 All means are the same 
#      H1 at least one mean is different

df_stats <- simple_df_agl[simple_df_agl$combination == "AGL1_D36E+HopM1",]
anova_pst <- aov(Pst_CFU ~ Pst_exp_date, data = df_stats)

{summary(anova_pst)
  summary.lm(anova_pst)}

# 4.2 TUKEY Posthoc to see which groups are significantly different 
TUKEY <- TukeyHSD(anova_pst)

# 4.3 Summarize TUKEY test results for saving & plotting

# 4.3.1 write a funktion to generate labels for the test results
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc results
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # order labels according to boxplot:
  Tukey.labels$Pst_exp_date = rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Pst_exp_date) , ]
  return(Tukey.labels)
}

# 4.3.2 Apply the function on the df
LABELS <- generate_label_df(TUKEY , "Pst_exp_date")
names(LABELS) <- c("Letters","Pst_strain")

# 4.3.3 safe TUKEY test results
sink(paste(pre, sep = "","ATTA_CFU_count_Agl-AGL1_D36E+HopM1.txt"))
print(LABELS)
sink()

###############################################################################

pst_raw <- read_xlsx("/home/melanie/working_directory/D36E_assays/ATTA_spinach/LDM20230524_CFU_count_with_DC3000.xlsx")

# 2. Select/ filter data ----
# 2.1 Bacteria strains as factors
pst_raw$Bacteria.strain <- as.factor(pst_raw$Bacterial_strain)
pst_raw$Combination <- as.factor(pst_raw$Combination)
#pst_raw$Date <- as.factor(pst_raw$dpi)
pst_raw$Date <- as.factor(
  ifelse(pst_raw$dpi == -2, "AminusTwo",
         ifelse(pst_raw$dpi == 0, "Bzero",
                ifelse(pst_raw$dpi == 3, "Cthree", pst_raw$dpi)))
)

# remove rows with NA (if extra rows are loaded)
# pst_raw <- pst_raw[1:210,]

# 2.2 log10 transformation of data for better presentation

#MIN <- min(pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`[pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`!=0]) # set lowest CFU counts as 0
pst_raw$CFU.log10 <- log10(pst_raw$`CFU/ml/1cm^2`+1)


# 2.3 make df containing only data I want to plot: Pst strain, CFU log10 transormed
Pst_strain <- as.factor(pst_raw$Bacteria.strain)
combination <- pst_raw$Combination
Pst_CFU <- pst_raw$CFU.log10
Pst_exp_date <- as.factor(pst_raw$Date)
ab <- pst_raw$`antibiotic plate`
#Pst_OD <- pst_raw$OD6000.infiltrated.0dpi
#Pst_comment <- pst_raw$Comment

simple_df <- data.frame(Pst_strain, combination, Pst_exp_date, Pst_CFU, ab)
simple_df <- subset(simple_df, Pst_strain %in% c("D36E", "AGL1_D36E", "DC3000", "AGL1_DC3000"))

simple_df_pst <- simple_df[simple_df$Pst_strain == "Pseudomonas",]
simple_df_agl <- simple_df[simple_df$Pst_strain == "AGL1",]
simple_df_two <- simple_df[simple_df$Pst_strain == "AGL1_and_Pseudomonas",]

# 2.4 Filter for group of interest
# simple_df_agl <- simple_df[simple_df$Pst_strain == "AGL1",]
# simple_df_agl <- simple_df_agl[simple_df_agl$Pst_CFU != 0, ]

co <- c("#5b5d5cff", "#44AA99", "#E69F00")  # Define colors
co <- c("#44AA99", "#E69F00")  # Define colors

p_agl <- ggplot(simple_df_pst, aes(x = as.factor(combination), y = Pst_CFU, fill = Pst_exp_date)) +
  geom_boxplot(outlier.colour = "red", alpha = 0.9, 
               color = "black",  # Black outlines, whiskers, and median
               fatten = 2) +  # Thicker median line
  labs(x = "Pst_exp_date", y = "log10 (CFU/cm2)") +
  geom_point(aes(color = Pst_exp_date), position = position_jitterdodge(0.1)) +  # Scatter colored by date
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_fill_manual(values = c("white", "white", "white")) +  # Fill colors for boxes
  scale_color_manual(values = co) +  # Colors for scatter points
  ylim(0, 9)

p_agl

ggsave(filename =  paste0("ATTA_CFU_count_Agl_boxplot.svg"), 
       plot = p_agl,
       device = "svg")

