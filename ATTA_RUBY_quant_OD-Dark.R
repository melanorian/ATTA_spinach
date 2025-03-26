# A) General basics ----
# 1.1 clean environment 
rm(list=ls())

# 1.2 generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")

# 2.1 open libraries
library(readxl)
library(ggplot2)
library(dplyr)
# library(tidyverse)
# library(plotly)
# library(tidyr)

# 2.2 set working directory
setwd("/home/melanie/working_directory/D36E_assays/ATTA_spinach/")

# 2.3. load file
df <- read_xlsx("LDM20221103_Dark_treatment_OD_repeat.xlsx")
df <- as.data.frame(df)

# 3. Boxplot with scatter

f_boxplot <- function(df, x_data, y_data) {
  
  # Create a boxplot with scatter overlay
  p <- ggplot(df, aes(x = as.factor(x_data), y = y_data, color = as.factor(OD))) +
    
    # Boxplot with black outline, no fill (so whiskers remain visible)
    geom_boxplot(outlier.shape = NA, color = "black", fill = NA, position = "identity") + 
    
    # Scatter points with jitter (larger dots for visibility)
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 4, alpha = 0.7) +
    
    theme_classic() +
    labs(x = "OD", y = "Percentage Leaf Area", color = "OD") +
    
    # Custom colors for scatter points and boxplots
    scale_color_manual(values = c("#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00")) +
    
    coord_cartesian(ylim = c(0, 30))
  
  return(p)
}

# Split data into "Dark" and "Light" subsets
df_list <- split(df, df$experiment)

# Apply f_boxplot() to each subset using lapply()
plots <- lapply(df_list, function(sub_df) {
  f_boxplot(sub_df, x_data = sub_df$experiment_OD_OD, y_data = sub_df$percentage_leaf_area)
})

# Assign plots to sperate variables
p_dark <- plots$Dark  # Plot for "Dark" experiment
p_light <- plots$Light # Plot for "Light" experiment

# Explore distribution - Os
# Example histogram of 'percentage_leaf_area' in your dataframe
df_no_control <-lapply(df_list, function(x) x[!grepl("OD_0$", x$experiment_OD_OD), ])

p_hist <- ggplot(df_no_control$Dark, aes(x = percentage_leaf_area)) +
  geom_histogram(binwidth = 0.5, fill = "#E69F00", color = "black", alpha = 1.0) +
 coord_cartesian(xlim = c(0, 28), ylim = c(0, 45)) +
  theme_minimal() +
  labs(title = "Histogram of Percentage Leaf Area", x = "Percentage Leaf Area", y = "Frequency")

# # Safe plots
ggsave(filename =  paste0("MM20250326_ATTA_RUBY_darkOD_his.svg"),
       plot = p_hist,
       device = "svg")

# Compare zero counts light/dark
zero_sum <- df_no_control$Light %>%
summarise(
count_zero = sum(percentage_leaf_area == 0),
count_between_0_and_1 = sum(percentage_leaf_area > 0 & percentage_leaf_area < 1),
total_count = n()
) %>%
mutate(
percent_zero = (count_zero / total_count) * 100,
percent_between_0_and_1 = (count_between_0_and_1 / total_count) * 100
)

# Total number of observations in each dataset
total_light <- 75
total_dark <- 77

# Number of values >1 (computed as remaining values)
count_above_1_light <- total_light - (25 + 26)  # Remaining values in OD_light
count_above_1_dark <- total_dark - (31 + 24)    # Remaining values in OD_dark

# Create contingency table
counts_matrix <- matrix(c(
  25, 26, count_above_1_light,  # OD_light
  31, 24, count_above_1_dark    # OD_dark
), nrow = 2, byrow = TRUE)

rownames(counts_matrix) <- c("OD_light", "OD_dark")
colnames(counts_matrix) <- c("Zero", "Between_0_and_1", ">1")

# Perform Chi-Square test
chi_test <- chisq.test(counts_matrix)

# Safe plots
# ggsave(filename =  paste0("ATTA_RUBY_lighOD_his.svg"),
#        plot = p_hist,
#        device = "svg")
# 
# ggsave(filename =  paste0("ATTA_RUBY_lightOD_boxplot.svg"), 
#        plot = p_light,
#        device = "svg")

# REMOVE zeros from dataset 
# Apply filtering to all dataframes in df_list
df_zerofree <- lapply(df_list, function(df) {
  df[df$percentage_leaf_area > 0 | grepl("OD_0$", df$experiment_OD_OD), ]
})

# Apply f_boxplot() to each subset using lapply()
plots_zerofree <- lapply(df_zerofree, function(sub_df) {
  f_boxplot(sub_df, x_data = sub_df$experiment_OD_OD, y_data = sub_df$percentage_leaf_area)
})

# Assign plots to sperate variables
p_dark <- plots_zerofree$Dark  # Plot for "Dark" experiment
p_light <- plots_zerofree$Light # Plot for "Light" experiment

ggsave(filename =  paste0("ATTA_RUBY_darkOD_noZero_box.svg"),
       plot = p_dark,
       device = "svg")

# 4. Statistical analysis across samples using ANOVA+TUKEY ----
# 4.1. ANOVA to analyse whether there are differences between the groups
#      H0 All means are the same 
#      H1 at least one mean is different

df_stats <- df_zerofree$Light
anova_pst <- aov(percentage_leaf_area ~ experiment_OD_OD, data = df_stats)

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
  Tukey.labels$experiment_OD_OD = rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$experiment_OD_OD) , ]
  return(Tukey.labels)
}

# 4.3.2 Apply the function on the df
LABELS <- generate_label_df(TUKEY , "experiment_OD_OD")
names(LABELS) <- c("Letters","Treatment")

# 4.3.3 safe TUKEY test results
sink(paste(pre, sep = "","ATTA_RUBY_OD_Light_noZero.txt"))
print(LABELS)
sink()

################################################################################
############################## Light vs Dark ###################################
################################################################################
df_zerofree <-df[df$percentage_leaf_area > 0 | grepl("OD_0$", df$experiment_OD_OD), ]

co <- c("#44AA99", "#E69F00") # define 2 colours 


# Grouped boxplot using ggplot
p_grouped <- ggplot(df_zerofree, aes(x = as.factor(OD), y = percentage_leaf_area, color = experiment, fill = experiment)) +
  geom_boxplot(outlier.colour = "white", alpha = 0, position = position_dodge(width = 0.75)) +  # Group boxplots side by side
  labs(x = "OD", y = "RUBY %") +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +  # Add jittered points
  theme_classic() +
  guides(x = guide_axis(angle = 45)) +  # Rotate x-axis labels for clarity
  scale_fill_manual(values = co) +  # Apply custom colors to fill
  scale_color_manual(values = co)   # Apply custom colors to color


ggsave(filename =  paste0("MM20250326_ATTA_RUBY_darkVSlight_boxplot.svg"),
       plot = p_grouped,
       device = "svg")

# 5. Statistical analysis that compares for each OD, Pst strain 0dpi and 2dpi (or days of interest)
# 5. Statistical analysis using a two-sided t-test ---- 

# 5.1 DEFINE VARIABLES OF INTEREST
OD <- unique(df_zerofree$OD) # extract the ODs from the orginal df
d1 <- "Light"  # Defines the reference dpi (usally 0dpi), epxeriment_date column as input  
d2 <- "Dark"  # Defines dpi of interest, usually 2dpi, epxeriment_date column as input

subset_df <- split(df_zerofree, df_zerofree$OD)

# df <- subset_df$'0'
# date1 <- d1
# date2 <- d2

tt <- function(df, OD, date1, date2) {
  a <- as.data.frame(df)
  gr1 <- a[(a$experiment == date1), ]
  gr2 <- a[(a$experiment == date2), ]
  b <- t.test(gr1$percentage_leaf_area, gr2$percentage_leaf_area, alternative = "two.sided", var.equal = TRUE)
  return(b$p.value)
}

result_list <- lapply(subset_df, tt, date1 = d1, date2 = d2)
df_t_res <- unlist(result_list)

# # 5.8 Update the significance column based on the p-value
# df_t_res[df_t_res$p_value > 0.05, "significance"] <- "n.s"
# df_t_res[df_t_res$p_value < 0.05, "significance"] <- "*"
# df_t_res[df_t_res$p_value < 0.01, "significance"] <- "**"
# df_t_res[df_t_res$p_value < 0.001, "significance"] <- "***"
# 
# # 5.9 safe the output as a .txt file
# sink(paste(pre, sep = "","_tt_summary_CUF_Pst_strains.txt"))
# print(df_t_res)
# sink()

################################################################################

# Example data (you can use your actual data frame)
df_light <- df_list$Light
df_dark <- df_list$Dark

# Step 2: Remove the prefix ("Light_" or "Dark_") from the "experiment_OD_OD" to match the groups
df_light <- df_light %>%
  mutate(experiment_OD_OD = gsub("^Light_", "", experiment_OD_OD))

df_dark <- df_dark %>%
  mutate(experiment_OD_OD = gsub("^Dark_", "", experiment_OD_OD))

# Step 1: Calculate mean per "experiment_OD_OD" group in the Light data
mean_light <- df_light %>%
  group_by(experiment_OD_OD) %>%
  summarise(mean_light = mean(percentage_leaf_area, na.rm = TRUE))

# Step 3: Merge the Light mean values with the Dark data based on "experiment_OD_OD"
df_dark_with_mean <- df_dark %>%
  left_join(mean_light, by = "experiment_OD_OD")

# Step 4: Calculate the percentage for each "Dark" value based on the corresponding "Light" mean
df_dark_with_mean <- df_dark_with_mean %>%
  mutate(percentage_dark = percentage_leaf_area-mean_light)  #(percentage_leaf_area / mean_light) * 100)

# # Step 5: Combine the "Light" and "Dark" data into a single dataframe
# df_combined <- bind_rows(df_light %>%
#                            mutate(experiment = "Light", percentage_leaf_area = percentage_leaf_area),
#                          df_dark_with_mean %>%
#                            mutate(experiment = "Dark", percentage_leaf_area = percentage_dark))

# Now, we can use your existing boxplot function
p_dvsl <- f_boxplot(df_dark_with_mean, df_dark_with_mean$experiment_OD_OD, df_dark_with_mean$percentage_dark)

# ggsave(filename =  paste0("ATTA_RUBY_darkVSlight_boxplot.svg"), 
#        plot = p_dvsl,
#        device = "svg")

################################################################################
################################### Model ######################################
################################################################################

df$Repeat <- rep("two", times = nrow(df))
#df$Condition <- gsub("_.*", "", df$experiment_OD_OD)

ggplot(df, aes(x = Repeat, y = percentage_leaf_area)) +
geom_point() + 
geom_smooth(method = "lm") +  # This shows a linear fit
theme_minimal()

# Example histogram of 'percentage_leaf_area' in your dataframe
p_hist <- ggplot(df, aes(x = percentage_leaf_area)) +
geom_histogram(binwidth = 0.1, fill = "#E69F00", color = "black", alpha = 0.7) +
theme_minimal() +
labs(title = "Histogram of Percentage Leaf Area", x = "Percentage Leaf Area", y = "Frequency")

# ggsave(filename =  paste0("ATTA_RUBY_darkVSlight_hist_bin01.svg"),
#        plot = p_hist,
#        device = "svg")

# Count how many zero values are present in the percentage_leaf_area column
sum(df$percentage_leaf_area == 0)  # Count of zeros
nrow(df)  # Total number of rows
sum(df$percentage_leaf_area == 0) / nrow(df) * 100



######################## Model expression Continouts ###########################
# A) Continous Model (Model Percentage Expression)
install.packages("tweedie")
library(tweedie)

# Fit a Tweedie model (power = 1.5 is a typical value for zero-inflated continuous data)
continuous_model <- glm(percentage_leaf_area ~ treatment + OD_infiltration + experiment_OD_OD,
                        data = df_combined, family = tweedie(link = "log", power = 1.5))

continuous_model <- glm(percentage_leaf_area ~ treatment + OD_infiltration + experiment_OD_OD,
                        data = df_combined, family = tweedie(p = 1.5))

# Model summary
summary(continuous_model)

################33

# Filter out zero values for continuous modeling
df_non_zero <- df[df$OD > 0, ]
df_gamma <- df_non_zero[df_non_zero$percentage_leaf_area  > 0, ]

# Fit a Gamma GLM with log link
gamma_model <- glm(percentage_leaf_area ~ experiment + OD,
                 data = df_gamma, family = Gamma(link = "log"))

# Model summary
summary(gamma_model)

# Example histogram of 'percentage_leaf_area' in your dataframe
p_hist_nZ <- ggplot(df_gamma, aes(x = percentage_leaf_area)) +
  geom_histogram(binwidth = 0.1, fill = "#E69F00", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histogram of Percentage Leaf Area", x = "Percentage Leaf Area", y = "Frequency")

ggsave(filename =  paste0("ATTA_RUBY_darkVSlight_noZero_hist_bin01.svg"),
       plot = p_hist_nZ,
       device = "svg")

## LM
# Fit a simple linear model for continuous percentage_leaf_area
lm_model <- lm(percentage_leaf_area ~ experiment + OD, data = df)

# Model summary
summary(lm_model)

# Count how many zero values are present in the percentage_leaf_area column
sum(df_non_zero$percentage_leaf_area == 0)  # Count of zeros
nrow(df_non_zero)  # Total number of rows
sum(df_non_zero$percentage_leaf_area == 0) / nrow(df_non_zero) * 100

########################## Model expression Logistic############################

# Create a new column 'Expression' based on whether percentage_leaf_area is zero or not
df$Expression <- ifelse(df$percentage_leaf_area > 0, 1, 0)

# Fit a logistic regression model to predict the binary outcome (Expression: TRUE/FALSE)
logistic_model <- glm(Expression ~ experiment + OD ,
                      data = df, family = "binomial")

# Model summary
summary(logistic_model)

# remove OD = 0 
df_no_zero_OD <- df[df$OD != 0, ]

# Fit a logistic regression model without OD = 0 values
logistic_model_no_zero_OD <- glm(Expression ~ experiment + OD,
                               data = df_no_zero_OD, family = "binomial")

# Model summary
summary(logistic_model_no_zero_OD)

# 4.3.3 safe TUKEY test results
sink(paste(pre, sep = "","ATTA_RUBY_OD_Dark_ZeroModel.txt"))
paste("Zero Modeling: OD 0 included")
summary(logistic_model)
paste("Zero Modeling: OD 0 excluded")
summary(logistic_model_no_zero_OD)
sink()


# Visualise? 
# Generate predicted values from the model
df_gamma$predicted_values <- predict(gamma_model, type = "response")

ggplot(df_gamma, aes(x = percentage_leaf_area, y = predicted_values)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Observed vs Predicted Values",
       x = "Observed Percentage of Leaf Area",
       y = "Predicted Percentage of Leaf Area") +
  theme_minimal()

# Predict values for a range of "experiment" while keeping "OD" constant
new_data <- data.frame(experiment = unique(df_gamma$experiment),
                       OD = mean(df_gamma$OD))

# Predict using the model
new_data$predicted_values <- predict(gamma_model, newdata = new_data, type = "response")

# Plot effect of experiment
ggplot(new_data, aes(x = experiment, y = predicted_values)) +
  geom_point(color = "green") +
  geom_line(color = "green") +
  labs(title = "Effect of Experiment on Predicted Leaf Area",
       x = "Experiment", y = "Predicted Percentage of Leaf Area") +
  theme_minimal()


# Generate predicted probabilities
df_no_zero_OD$predicted_probabilities <- predict(logistic_model_no_zero_OD, type = "response")


ggplot(df_no_zero_OD, aes(x = Expression, y = predicted_probabilities)) +
  geom_jitter(width = 0.1, height = 0.1, color = "blue", alpha = 0.5) + 
  geom_smooth(method = "loess", color = "red") +
  labs(title = "Observed vs Predicted Probabilities",
       x = "Observed Expression (0 or 1)", 
       y = "Predicted Probability of Expression") +
  theme_minimal()
