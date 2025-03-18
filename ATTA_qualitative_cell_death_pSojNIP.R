# A) General basics ----
# 1. clean environment 
rm(list=ls())

# generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")

# 2. open libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readxl)
library(plotly)
library(tidyr)

# 3. set working directory
setwd("/home/melanie/working_directory/D36E_assays/disease_data/")

# 4. load file
df <- read_xlsx("LDM20221129_RUBY_psojNIP_Experiment.xlsx")
# df <- as.data.frame(df)
# rownames(df) <- df$category
# df$category <- NULL

# B) get data(frame) into desired shape ----
# 5. rename columns and select only relevant columns
df1 <- data.frame("category" = df$category,
                 "nr_score_0" = df$score_1,
                 "nr_score_1" = df$score_2,
                 "nr_score_2" = df$score_3,
                 "nr_score_3" = df$score_4,
                 "nr_score_4" = df$score_5,
                 "nr_score_5" = df$score_6
)

# rownames(df1) <- df$category
# df1$category <- NULL

# Separate the 'category' column into four parts: OD, Number, Leaf, and Treatment
df1 <- df1 %>%
  separate(category, into = c("OD", "Number", "Leaf", "Treatment"), sep = "_", extra = "merge")

df1$OD <- NULL

##############
##############
# Convert the wide format into long format
df1_long <- df1 %>%
  pivot_longer(cols = starts_with("nr_score_"), 
               names_to = "Score", 
               values_to = "Count")

# Clean up the "Score" column (remove "nr_score_" prefix)
df1_long$Score <- factor(gsub("nr_score_", "Score ", df1_long$Score), 
                         levels = c("Score 0", "Score 1", "Score 2", "Score 3", "Score 4", "Score 5"))


custom_colors <- c(
  "Score 0" = "#fffbf1",  # Orange
  "Score 1" = "#ffeabc",  # Blue
  "Score 2" = "#ffbe2d",  # Green
  "Score 3" = "#E69F00",  # Yellow
  "Score 4" = "#dd7e00ff",  # Dark Blue
  "Score 5" = "#D55E00"    # Red
)

#c("#fff2d6ff", "#fff2d6ff", "#ffeabcff", "#ffd374ff", "#ffbe2dff","#e69f00ff", "#dd7e00ff", "#d55e00ff", "#9f4400ff")

# Create stacked bar graph with custom colors
ggplot(df1_long, aes(x = Treatment, y = Count, fill = Score)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  theme_classic() +
  facet_grid(Number ~ Leaf) +  # Creates separate plots for each (Number, Leaf) combination
  labs(x = "Treatment", y = "Count", fill = "Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_fill_manual(values = custom_colors)  # Apply custom colors


ggsave(paste0("/home/melanie/working_directory/D36E_assays/ATTA_spinach/", pre, "_scorePsojNIP.svg"),width = 12, height = 10)

############################################################################# 
############################# Statistics ####################################
#############################################################################

# E) Prepare data for Stats (partly redundand due to recycling of script) ----

# df with relevant information. Go back to original df 
### calculate the proportions in % of cell-death positive/negative infiltration spots
df2 <- data.frame("category" = df$category,
                  "Effector_line" = df$`Effector line`, 
                  "infiltration_spots" = df$`n (infiltration spots)`,
                  "disease_pos" = df$`Positive = visual signs of cell death`,
                  "disease_neg" = df$`Negative = no visual signs of cell death`)
                  

# drop first row that only contains numbers of disease scores
df2 <- df2[!(rownames(df2) == '1'),]                 

# write function to sum up all the +/- infiltration spots (to account for repeating control samples)

sum_d <-  lapply(e_lines, function(x){ # sum up of cell-death POSITIVE
  a <- df2[(df2$Effector_line == x),] # subset datset by effector line
  b <- sum(a$disease_pos) # sum infiltration spots WITH cell death
  return(b)
}) %>% unlist()

sum_n <- lapply(e_lines, function(x){ # sum up of cell-death NEGATIVE
  a <- df2[(df2$Effector_line == x),] # subset datset by effector line
  b <- sum(a$disease_neg) # sum infiltration spots NO cell death
  return(b)
}) %>% unlist()


# store in new dataframe
simple_df_d <- as.data.frame(cbind(unlist(e_lines), as.numeric(sum_d))) # df cell death
simple_df_n <- as.data.frame(cbind(unlist(e_lines), as.numeric(sum_n))) # df no cell death


# rename columns of simple_df in a sensible way
colnames(simple_df_d)[c(1, 2)] = c("Effector", "count_death")
colnames(simple_df_n)[c(1, 2)] = c("Effector", "count_death")

# add column with "condition" which symptom positive or negative 
# this is necessary to get data into the right shape for plotting
simple_df_d["condition"] <- as.vector(rep("death", times = length(e_lines)))
simple_df_n["condition"] <- as.vector(rep("no death", times = length(e_lines)))

# join the data frame for plotting
simple_df_2 <- full_join(simple_df_d, simple_df_n)


# function to sum up 
nr_spots <- lapply(e_lines, function(x){
  line <- simple_df_2[simple_df_2$Effector == x,]
  s <- sum(as.numeric(line$count_death))
  return(s)
}) 

# add nr all infiltration spots to df
simple_df_2["nr_spots"] <- unlist(nr_spots) 

# calculate percentag
percentage <- (as.numeric(simple_df_2$count_death)/as.numeric(simple_df_2$nr_spots))*100
simple_df_2["Percentage"] <- percentage


# F) # Statistical analysis of results. Chi-square test cannot be performed due to too few
# samples (<5) in some groups --> Fisher's exact test ----

# bring datafram in a shape that's suitable for Fisher's exact test 
# re-name df so that info on cell-death-pos/neg can be merged based on effectors


simple_df_d$condition <- NULL
colnames(simple_df_d) <- c("Effector", "cell-death-positive")


simple_df_n$condition <- NULL
colnames(simple_df_n) <- c("Effector", "cell-death-negative")


# construct df for stat analysis
df_stat <- merge(simple_df_d, simple_df_n)

# Effectors as colname
rownames(df_stat) <- df_stat$Effector
df_stat$Effector <- NULL

# factors as numeric for applying statistics (oddly complicated but well...)

i <- c(1, 2)
df_stat[ , i] <- apply(df_stat[ , i], 2,          
                       function(x) as.numeric(as.character(x))) # first convert to character, important to retain the values


effect <- unique(rownames(df_stat))

# create a frequency table
#count(df_stat, df_stat$`cell-death-positive`)

#typeof(effect)
#df_stat <- as.factor(colnames(df_stat))

# expected frequencies per effector
#chisq.test(df_stat)$expected # --> results in error due to low expected frequencies

F_test <- sapply(effect, function(x){
  d <- df_stat[x,]                                #subset df to extract single eff row
  c <- df_stat[(rownames(df_stat) == "AAAD36E"),] #subset df to extract control row
  m <- rbind(c ,d)                                #join both rows to have a small df
  f <- fisher.test(m)                             #apply test
  p <- f$p.value
  return(p)}) %>% cbind() %>% data.frame()

colnames(F_test) <- "p-value"

# creat new column containing values of row names so that we have an identical
# column in both data sets that we can merge on 
F_test$Effectors <- row.names(F_test)
df_stat$Effectors <- row.names(df_stat)

# merge datafram

m <- merge(df_stat, F_test)

# add column with % of disease-symptom positive infiltration spots
perc_pos <- subset(simple_df_2, !grepl("no death", condition))
perc_pos <- arrange(perc_pos, Effector )

m$perc_pos <- perc_pos$Percentage

# safe results in a txt file
write.table(m, file = paste(pre, sep = "","_stats_cell_death_2_3_4_positive.txt"))

##############################################################################
#########################  export for summary heat map #######################
##############################################################################

# list all unique strains
e_lines

# only keep % of cell death positive incidence
d2 <- simple_df[!(simple_df$condition == "no death"),]

# drop all the columns not needed
d2 <- d2[c("Effector","Percentage")]

# Rename D36E, DC3000 to appear in the same order outputs in other analysis
d2[d2 == "AAAD36E"] <- "D36E"
d2[d2 == "ZZZDC3000"] <- "DC3000"

# export 
write.csv(d2, "./2_Cell_death_positive_percentage.csv")
