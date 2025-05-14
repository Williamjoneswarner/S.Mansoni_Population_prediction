#Prevalence Prediction from KK intentisty of infection using linear regression. 

#Utilisation of a model trained on KK1 from a single stool to predict probability of detecting a positive
#at this EPG based on a village population also tested with 4KK from 4 stools. Detecting Schistosoma Mansoni

########Data Manipulation########

df <- read.csv("SM_Data_4kk_n327.csv")
dim(df)
head(df)

#Selecting only Columns required
selected_df <- df[, c(2:6, 8, 10)]
head(selected_df)

#Multiplying to get EPG - not needed but helps visualise the data. 
selected_df[1:4] <- selected_df[1:4]*20
View(selected_df)
# Load the tidyverse package
library(tidyverse)

# Use pivot_longer to reshape the data - treating each KK1 for each individual as a seporate event. 
#maximising number of observations and the variability possible for each KK4 event. 
long_df <- selected_df %>%
  pivot_longer(
    cols = -c(MPEG, Prob_pos, KK_POS),  # Exclude both KK_4_EPG and KK_POS
    names_to = "Combination",
    values_to = "Value"
  )

# View the first few rows of the long format dataframe
head(long_df)

library(dplyr)
selected_df <- select(long_df, -Combination)

head(selected_df)

selected_df <- selected_df %>%
  rename(KK_1_EPG = Value)

selected_df <- selected_df %>%
  rename(KK_4_EPG = MPEG)

selected_df <- selected_df[, c("KK_4_EPG", "KK_1_EPG", "KK_POS", "Prob_pos")]

write.csv(selected_df, "KK_Count_long_4KK.csv")

head(selected_df)
dim(selected_df)

combined_df <- selected_df

# Calculate Spearman correlation and p-value
cor_test <- cor.test(combined_df$KK_1_EPG, combined_df$KK_4_EPG, method = "spearman")

# Extract the Spearman correlation coefficient
correlation <- cor_test$estimate

library(ggplot2)

# Create ggplot
ggplot(combined_df, aes(x = KK_1_EPG, y = KK_4_EPG)) +
  geom_point(color = "blue") +  # Plot the points
  geom_smooth(method = "lm", se = TRUE, color = "red") +  # Add linear model line
  labs(
    title = "Intensity of Infection Compared (2 vs 6 KK)",
    x = "Intensity (1KK EPG)",
    y = "Intensity (4KK EPG)"
  ) +
  annotate("text", x = 1000, y = 10, 
           label = paste("Spearman Correlation: ", round(correlation, 2)), 
           color = "black", size = 5)

##creating table for analysis - grouping by positive KK's by KK1 EPG. 

# Convert table to dataframe
library(dplyr)

# Create df_table using dplyr instead of table()
df_table <- combined_df %>%
  count(KK_POS, KK_1_EPG, name = "Count") %>%  # Count occurrences
  arrange(KK_1_EPG, KK_POS)  # Optional: Arrange for better readability

# Print the dataframe
print(df_table)

library(dplyr)

# Compute weighted sum and total count by KK_1_EPG - probabiltiy of positive KK by each EPG 
# Create the summarized table with counts for each KK_POS as columns
df_table_sum <- df_table %>%
  mutate(weighted_value = KK_POS * Count) %>%  # Calculate weighted value
  group_by(KK_1_EPG) %>%  # Group by KK_1_EPG
  summarise(
    Sum_Weighted_Value = sum(weighted_value),  # Sum of weighted values
    Total_Count = sum(Count) * 4,  # Sum of counts for each KK_1_EPG
    `0` = sum(Count[KK_POS == 0]),  # Count for KK_POS == 0
    `1` = sum(Count[KK_POS == 1]),  # Count for KK_POS == 1
    `2` = sum(Count[KK_POS == 2]),  # Count for KK_POS == 2
    `3` = sum(Count[KK_POS == 3]),  # Count for KK_POS == 2
    `4` = sum(Count[KK_POS == 4])   # Count for KK_POS == 3
  ) %>%
  ungroup()  # Ensure no grouping remains

# Print the result
print(df_table_sum)

# Assuming df_table_sum already exists as shown in the example - calculate proportion of positive 
#KK's at each reported KK1 EPG

df_table_sum <- df_table_sum %>%
  mutate(Ratio_Weighted_Value_Per_Count = Sum_Weighted_Value / Total_Count)

# Print the updated dataframe
print(df_table_sum)

# Create the plot
ggplot(df_table_sum, aes(x = KK_1_EPG, y = Ratio_Weighted_Value_Per_Count)) +
  geom_line() +  # Line plot
  geom_point() +  # Add points to the plot
  labs(
    title = "Ratio of Positive/Negative KK slides at Given EPG",
    x = "Intesity of Infection (2KK)",
    y = "Ratio of Positive KK"
  )

# Perform the join and merge the probability values back with 
combined_df <- combined_df %>%
  left_join(df_table_sum %>%
              select(KK_1_EPG, Ratio_Weighted_Value_Per_Count), 
            by = "KK_1_EPG")

# View the updated combined_df
head(combined_df)

### building the model####
#logit transform the probabity of identiying the kato katz at a speicfic EPG
# Ensure the values are between 0 and 1, if necessary
# Adjust values that are exactly 0 or 1 to avoid logit issues

combined_df$Ratio_Weighted_Value_Per_Count <- as.numeric(combined_df$Ratio_Weighted_Value_Per_Count)

combined_df <- combined_df %>%
  mutate(
    Ratio_adjusted = ifelse(Ratio_Weighted_Value_Per_Count == 1, 0.9999999,
                            Ratio_Weighted_Value_Per_Count),
    Ratio_logit = log(Ratio_adjusted / (1 - Ratio_adjusted))
  )

#Check its worked. 
max(combined_df$Ratio_adjusted)
min(combined_df$Ratio_adjusted)

model_logit <- lm(Ratio_logit ~ KK_1_EPG, data = combined_df)
model_logit_lm <- lm(Ratio_logit ~ sqrt(KK_1_EPG), data = combined_df)


anova(model_logit_lm, model_logit)

# Step 2: Get predicted logit_base values for basic model
combined_df$predicted_logit_base <- predict(model_logit, newdata = combined_df)

# Reverse logit transformation to get predicted values back on the original scale
combined_df$predicted_base_original <- 1 / (1 + exp(-combined_df$predicted_logit_base))

ggplot(combined_df, aes(x = KK_1_EPG, y = Ratio_Weighted_Value_Per_Count)) +
  geom_point(alpha = 0.6) +  # Observed points
  geom_line(aes(x = KK_1_EPG, y = predicted_base_original), color = "blue") +  # Predicted line
  labs(
    title = "Observed vs Predicted Values (Blue = sqrt(KK_1_EPG)",
    x = "KK_1_EPG",
    y = "Probability of a Positive Kato-Ktaz"
  )

# Step 2: Get predicted logit values for sqrt() model
combined_df$predicted_logit <- predict(model_logit_lm, newdata = combined_df)

# Reverse logit transformation to get predicted values back on the original scale
combined_df$predicted_original <- 1 / (1 + exp(-combined_df$predicted_logit))

ggplot(combined_df, aes(x = KK_1_EPG, y = Ratio_Weighted_Value_Per_Count)) +
  geom_point(alpha = 0.6) +  # Observed points
  geom_line(aes(x = KK_1_EPG, y = predicted_original), color = "blue") +  # Predicted line
  labs(
    title = "Observed vs Predicted Values (Blue = sqrt(KK_1_EPG)",
    x = "KK_1_EPG",
    y = "Probability of a Positive Kato-Ktaz"
  )

#####Additional Models exploring GLM and Generalized Additive Model######
View(df_table_sum)

#this uses the values obtained in the summary table 'df_table_sum'
model_glm <- glm(cbind(Sum_Weighted_Value, Total_Count - Sum_Weighted_Value) ~ sqrt(KK_1_EPG),
                 family = binomial(link = "logit"), data = df_table_sum)

df_table_sum$predicted_prob <- predict(model_glm, type = "response")

#now join with the original dataset with individual level data
library(dplyr)

combined_df <- combined_df %>%
  left_join(df_table_sum %>% select(KK_1_EPG, Sum_Weighted_Value, Total_Count), by = "KK_1_EPG")

combined_df <- combined_df %>%
  select(-Sum_Weighted_Value, -Total_Count)

combined_df$predicted_prob_glm <- predict(model_glm, newdata = combined_df, type = "response")

library(ggplot2)

ggplot(df_table_sum, aes(x = sqrt(KK_1_EPG))) +
  geom_point(aes(y = Sum_Weighted_Value / Total_Count), color = "black", size = 1) +  # Observed proportion
  geom_line(aes(y = predicted_prob), color = "blue", size = 0.7) +  # Predicted probability
  labs(
    x = "sqrt(KK_1_EPG)",
    y = "Proportion Positive",
    title = "Observed vs Predicted Proportion Positive"
  )

#GAM model
library(mgcv)

model_gam <- gam(Ratio_logit ~ s(KK_1_EPG), data = combined_df)

combined_df$predicted_logit_gam <- predict(model_gam)

combined_df$predicted_prob_gam <- 1 / (1 + exp(-combined_df$predicted_logit_gam))

#plot all models on single plot with raw data. 
ggplot(combined_df, aes(x = KK_1_EPG)) +
  geom_point(aes(y = Ratio_adjusted, color = "Observed"), size = 1) +
  geom_line(aes(y = predicted_prob_gam, color = "GAM Prediction"), size = 0.7) +
  geom_line(aes(y = predicted_original, color = "Logit LM Prediction (sqrt)"), size = 0.7) +
  geom_line(aes(y = predicted_prob_glm, color = "GLM Prediction"), size = 0.7) +
  geom_line(aes(y = predicted_base_original, color = "Logit LM Prediction (linear)"), size = 0.7) +
  scale_color_manual(
    name = "Legend",
    values = c(
      "Observed" = "black",
      "GAM Prediction" = "blue",
      "Logit LM Prediction (sqrt)" = "red",
      "GLM Prediction" = "orange",
      "Logit LM Prediction (linear)" = "green"
    )
  ) +
  labs(
    x = "KK_1_EPG",
    y = "Proportion Positive",
    title = "Observed vs Predicted Proportion Positive (All Models)"
  )

#model_logit
#model_logit_lm
#model_glm
#model_gam

######Testing model predictions######

#Reverse the probability of being detected to calculate the probability of missed diagnosis
#if this was detected how many were missed? 
combined_df$ProbPosScor_base <- 1/combined_df$predicted_base_original
combined_df$ProbPosScor <- 1/combined_df$predicted_original
combined_df$ProbPosScor_gam <- 1/combined_df$predicted_prob_gam
combined_df$ProbPosScor_glm <- 1/combined_df$predicted_prob_glm

#only calculate on the positive individuals.
#this is the number of probable positives in the population based off the individual level EPG divided by the numbr of individuals in the population
pred_prev_base <- sum(combined_df$ProbPosScor_base[combined_df$KK_1_EPG != 0]) / nrow(combined_df)
pred_prev <- sum(combined_df$ProbPosScor[combined_df$KK_1_EPG != 0]) / nrow(combined_df)
pred_prev_gam <- sum(combined_df$ProbPosScor_gam[combined_df$KK_1_EPG != 0]) / nrow(combined_df)
pred_prev_glm <- sum(combined_df$ProbPosScor_glm[combined_df$KK_1_EPG != 0]) / nrow(combined_df)

#View the predicted prevalence of positive individuals
pred_prev_base
#0.7815509
pred_prev
#0.7018444
pred_prev_gam
#0.8471201
pred_prev_glm
#0.7385028

#calculate the prevalence as determien by KK4
KK4_prev <- sum(combined_df$KK_4_EPG != 0) / nrow(combined_df)
#0.7607362

#####Test and Train of model same visits 1:10 split CV####
#this gives us the ability to explore how each model works at different prevalence levels and vague sense of accuracy on unknown data
names(combined_df)

set.seed(123)  # For reproducibility

run_model_multiple_times <- function(data, n_iter = 50) {
  results <- data.frame(
    iteration = integer(), 
    KK1_mean = numeric(),
    KK1_median = numeric(),
    KK1_max = numeric(),
    KK1_pos_count = integer(),
    light_count = integer(),
    moderate_count = integer(),
    heavy_count = integer(),
    pred_prev_logit = numeric(),
    pred_prev_base = numeric(),
    pred_prev_glm = numeric(),
    pred_prev_gam = numeric(),
    KK4_prev = numeric(),
    KK2_prev = numeric()
  )  
  
  for (i in 1:n_iter) {
    # Random split
    train_indices <- sample(seq_len(nrow(data)), size = 0.9 * nrow(data))
    train_data <- data[train_indices, ]
    test_data <- data[-train_indices, ]
    
    # Compute summary stats from test set
    KK1_mean <- mean(test_data$KK_1_EPG, na.rm = TRUE)
    KK1_median <- median(test_data$KK_1_EPG, na.rm = TRUE)
    KK1_max <- max(test_data$KK_1_EPG, na.rm = TRUE)
    KK1_pos_count <- sum(test_data$KK_1_EPG > 0, na.rm = TRUE)
    
    # Intensity categories
    light_count <- sum(test_data$KK_1_EPG >= 1 & test_data$KK_1_EPG <= 99, na.rm = TRUE)
    moderate_count <- sum(test_data$KK_1_EPG >= 100 & test_data$KK_1_EPG <= 399, na.rm = TRUE)
    heavy_count <- sum(test_data$KK_1_EPG >= 400, na.rm = TRUE)
    
    # Compute prevalence
    KK4_prev <- sum(test_data$KK_4_EPG != 0) / nrow(test_data)  
    KK2_prev <- sum(test_data$KK_1_EPG != 0) / nrow(test_data)  
    
    # Transform and fit models
    train_data <- train_data %>%
      mutate(
        Ratio_adjusted = pmin(pmax(Ratio_Weighted_Value_Per_Count, 1e-10), 1 - 0.000001),
        Ratio_logit = log(Ratio_adjusted / (1 - Ratio_adjusted))
      )
    
    test_data <- test_data %>%
      mutate(Ratio_adjusted = pmin(pmax(Ratio_Weighted_Value_Per_Count, 1e-10), 1 - 0.000001))
    
    model_logit_lm <- lm(Ratio_logit ~ sqrt(KK_1_EPG), data = train_data)
    model_base <- lm(Ratio_logit ~ KK_1_EPG, data = train_data)
    model_glm <- glm(Ratio_adjusted ~ sqrt(KK_1_EPG), data = train_data, family = binomial)
    model_gam <- gam(Ratio_adjusted ~ s(sqrt(KK_1_EPG)), data = train_data)
    
    # Predict
    test_data$predicted_logit <- predict(model_logit_lm, newdata = test_data)
    test_data$predicted_original_logit <- 1 / (1 + exp(-test_data$predicted_logit))
    
    test_data$predicted_base <- predict(model_base, newdata = test_data)
    test_data$predicted_original_base <- 1 / (1 + exp(-test_data$predicted_base))
    
    test_data$predicted_glm <- predict(model_glm, newdata = test_data, type = "response")
    test_data$predicted_gam <- predict(model_gam, newdata = test_data, type = "response")
    
    # Prevalence prediction
    pred_prev_logit <- sum(1 / test_data$predicted_original_logit[test_data$KK_1_EPG != 0]) / nrow(test_data)
    pred_prev_base <- sum(1 / test_data$predicted_original_base[test_data$KK_1_EPG != 0]) / nrow(test_data)
    pred_prev_glm <- sum(1 / test_data$predicted_glm[test_data$KK_1_EPG != 0]) / nrow(test_data)
    pred_prev_gam <- sum(1 / test_data$predicted_gam[test_data$KK_1_EPG != 0]) / nrow(test_data)
    
    # Store results
    results <- rbind(results, data.frame(
      iteration = i,
      KK1_mean = KK1_mean,
      KK1_median = KK1_median,
      KK1_max = KK1_max,
      KK1_pos_count = KK1_pos_count,
      light_count = light_count,
      moderate_count = moderate_count,
      heavy_count = heavy_count,
      pred_prev_logit = pred_prev_logit,
      pred_prev_base = pred_prev_base,
      pred_prev_glm = pred_prev_glm,
      pred_prev_gam = pred_prev_gam,
      KK4_prev = KK4_prev,
      KK2_prev = KK2_prev
    ))
  }
  
  return(results)
}

# Run the function
results_v1 <- run_model_multiple_times(combined_df)

View(results_v1)

final_results <- results_v1

# Define columns to summarise
cols_to_summarise <- c("pred_prev_logit", "pred_prev_base", "pred_prev_glm", 
                       "pred_prev_gam", "KK4_prev", "KK2_prev")

# Compute summary statistics
summary_results_iter <- data.frame(
  Variable = cols_to_summarise,
  Mean = sapply(results_v1[cols_to_summarise], mean),
  Lower_2.5_CI = sapply(results_v1[cols_to_summarise], function(x) quantile(x, 0.025)),
  Upper_97.5_CI = sapply(results_v1[cols_to_summarise], function(x) quantile(x, 0.975))
)

print(summary_results_iter)

write.csv(results_v1, "Multi_model_output.csv")

#now to visualise the results

library(tidyverse)

# Reshape data from wide to long format
results_long <- results_v1 %>%
  pivot_longer(cols = starts_with("pred_prev"), 
               names_to = "model", 
               values_to = "pred_prev") %>%
  mutate(model = recode(model,
                        pred_prev_logit = "Logit (sqrt)",
                        pred_prev_base = "Base (linear)",
                        pred_prev_glm = "GLM",
                        pred_prev_gam = "GAM"))

#plotting results of the train and hold out for the 4 models together. 
ggplot(results_long, aes(x = pred_prev, y = KK4_prev, color = model)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +  # Add linear model line per group
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y reference line
  labs(
    title = "Predicted vs. Observed KK4 Prevalence by Model",
    x = "Predicted Prevalence",
    y = "Observed KK4 Prevalence",
    color = "Model Type"
  ) +
  xlim(0.5, 1) +
  ylim(0.5, 1)

#now individually
ggplot(results_long, aes(x = pred_prev, y = KK4_prev, color = model)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Predicted vs. Observed KK4 Prevalence by Model (Faceted)",
    x = "Predicted Prevalence",
    y = "Observed KK4 Prevalence",
    color = "Model Type"
  ) +
  xlim(0.5, 1) +
  ylim(0.5, 1) +
  facet_wrap(~ model, ncol = 2)

# Step 1: Reshape and calculate pred_dif 
# Add difference from 'true' prevalence and predcited prevalence
results_long <- results_long %>%
  mutate(pred_dif = pred_prev - KK4_prev)


# Step 2: Plot pred_dif vs KK2_prev
ggplot(results_long, aes(x = KK2_prev, y = pred_dif, color = model)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ model, ncol = 2) +
  labs(
    title = "Prediction Difference (Predicted - KK4) vs KK2 Prevalence",
    x = "KK2 Prevalence",
    y = "Prediction Difference"
  )

#plot on same graph
ggplot(results_long, aes(x = KK2_prev, y = pred_dif, color = model)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Prediction Difference vs KK2 Prevalence (All Models)",
    x = "KK2 Prevalence",
    y = "Prediction Difference",
    color = "Model"
  )

###### Testing this model's approach on two new villages#####

#manipulating the new dataset to fit the models inputs. 
df <- read.csv("SM_Data.csv")

dim(df)

names(df)

df_2 <- df[df$Location == 2, ]
df_3 <- df[df$Location == 1, ]

# For df_2 2_kk
selected_df_2 <- df_2[, c(1:3, 6:8)]
selected_df_2[2:3] <- selected_df_2[2:3]*20
selected_df_2 <- selected_df_2[complete.cases(selected_df_2[, 1:3]), ]
selected_df_2$Prob_pos <- selected_df_2$KK_POS/selected_df_2$KK

# For df_3 3_kk
selected_df_3 <- df_3[, c(1:4, 6:8)]
selected_df_3[2:4] <- selected_df_3[2:4]*20
selected_df_3 <- selected_df_3[complete.cases(selected_df_3[, 1:5]), ]
selected_df_3$Prob_pos <- selected_df_3$KK_POS/selected_df_3$KK

library(tidyr)

# Use pivot_longer to reshape the data
long_df_2 <- selected_df_2 %>%
  pivot_longer(
    cols = -c(MPEG, Prob_pos, KK_POS, age, KK),  # Exclude both KK_4_EPG and KK_POS
    names_to = "Combination",
    values_to = "Value"
  )

# View the first few rows of the long format dataframe
head(long_df_2)

library(dplyr)
selected_df_2 <- select(long_df_2, -Combination)

head(selected_df_2)

selected_df_2 <- selected_df_2 %>%
  rename(KK_1_EPG = Value)

selected_df_2 <- selected_df_2 %>%
  rename(KK_2_EPG = MPEG)

selected_df_2 <- selected_df_2[, c("KK_2_EPG", "KK_1_EPG", "KK_POS", "Prob_pos")]

write.csv(selected_df_2, "KK_Count_long_2KK.csv")

# Use pivot_longer to reshape the data
long_df_3 <- selected_df_3 %>%
  pivot_longer(
    cols = -c(MPEG, Prob_pos, KK_POS, age, KK),  # Exclude both KK_4_EPG and KK_POS
    names_to = "Combination",
    values_to = "Value"
  )

# View the first few rows of the long format dataframe
head(long_df_3)

library(dplyr)
selected_df_3 <- select(long_df_3, -Combination)

head(selected_df_3)

selected_df_3 <- selected_df_3 %>%
  rename(KK_1_EPG = Value)

selected_df_3 <- selected_df_3 %>%
  rename(KK_3_EPG = MPEG)

selected_df_3 <- selected_df_3[, c("KK_3_EPG", "KK_1_EPG", "KK_POS", "Prob_pos")]

write.csv(selected_df_3, "KK_Count_long_3KK.csv")

#claculating the prevalence of the dataset if 1KK used and if 2KK/3KK used
KK_1_EPG_prev_value_2 <- sum(selected_df_2$KK_1_EPG > 0) / sum(selected_df_2$KK_1_EPG >= 0)
KK_2_EPG_prev_value_2 <- sum(selected_df_2$KK_2_EPG > 0) / sum(selected_df_2$KK_2_EPG >= 0)

KK_1_EPG_prev_value_2
KK_2_EPG_prev_value_2

KK_1_EPG_prev_value_3 <- sum(selected_df_3$KK_1_EPG > 0) / sum(selected_df_3$KK_1_EPG >= 0)
KK_3_EPG_prev_value_3 <- sum(selected_df_3$KK_3_EPG > 0) / sum(selected_df_3$KK_3_EPG >= 0)

KK_1_EPG_prev_value_3
KK_3_EPG_prev_value_3

combined_df_2 <- selected_df_2
combined_df_3 <- selected_df_3

#model_logit
#model_logit_lm
#model_glm
#model_gam

####Applying the models to the new dataset's####

ci_preds <- predict(model_logit_lm, newdata = combined_df_2, interval = "confidence")
# Back-transform CI predictions
ci_preds_original <- as.data.frame(ci_preds)
ci_preds_original <- ci_preds_original %>%
  mutate(
    fit = plogis(fit),
    lwr = plogis(lwr),
    upr = plogis(upr)
  )

#predicted and reverse logit
combined_df_2$predicted_logit <- predict(model_logit_lm, newdata = combined_df_2)
combined_df_2$predicted_original_logit <- 1 / (1 + exp(-combined_df_2$predicted_logit))

combined_df_2$predicted_base <- predict(model_logit, newdata = combined_df_2)
combined_df_2$predicted_original_base <- 1 / (1 + exp(-combined_df_2$predicted_base))

combined_df_2$predicted_glm <- predict(model_glm, newdata = combined_df_2, type = "response")

combined_df_2$predicted_gam <- predict(model_gam, newdata = combined_df_2, type = "response")
combined_df_2$predicted_original_gam <- 1 / (1 + exp(-combined_df_2$predicted_gam))

# Prevalence prediction
pred_prev_logit <- sum(1 / combined_df_2$predicted_original_logit[combined_df_2$KK_1_EPG != 0]) / nrow(combined_df_2)
pred_prev_base <- sum(1 / combined_df_2$predicted_original_base[combined_df_2$KK_1_EPG != 0]) / nrow(combined_df_2)
pred_prev_glm <- sum(1 / combined_df_2$predicted_glm[combined_df_2$KK_1_EPG != 0]) / nrow(combined_df_2)
pred_prev_gam <- sum(1 / combined_df_2$predicted_original_gam[combined_df_2$KK_1_EPG != 0]) / nrow(combined_df_2)

pred_prev_logit
pred_prev_base
pred_prev_glm
pred_prev_gam
#######

# ---- Load required packages ----
library(dplyr)
library(tibble)

# ---- Function to extract prediction summary ----
get_prediction_summary <- function(df, fit_col, lwr_col, upr_col) {
  data.frame(
    Prevalence = sum(1 / df[[fit_col]][df$KK_1_EPG > 0])/ nrow(df),
    CI_low = sum(1 / df[[lwr_col]][df$KK_1_EPG > 0]) / nrow(df),
    CI_high = sum(1 / df[[upr_col]][df$KK_1_EPG > 0]) / nrow(df)
  )
}



#sum(1 / combined_df_2$predicted_original_logit[combined_df_2$KK_1_EPG != 0]) / nrow(combined_df_2)

# ---- Create a new data frame with only the KK_1_EPG column ----
model_outputs_2 <- combined_df_2 %>%
  select(KK_1_EPG)

# ---- View the new data frame ----
head(model_outputs_2)

# ---- Predict with confidence intervals for each model ----
logit_lm_ci <- predict(model_logit_lm, newdata = model_outputs_2, interval = "confidence")

base_lm_ci  <- predict(model_logit, newdata = model_outputs_2, interval = "confidence")

# For GLM: type = "link" to get logit scale (like lm), then CI
glm_link_ci <- predict(model_glm, newdata = model_outputs_2, type = "link", se.fit = TRUE)
glm_ci_df <- data.frame(
  fit = glm_link_ci$fit,
  lwr = glm_link_ci$fit - 1.96 * glm_link_ci$se.fit,
  upr = glm_link_ci$fit + 1.96 * glm_link_ci$se.fit
)

# For GAM (same as GLM)
gam_link_ci <- predict(model_gam, newdata = model_outputs_2, type = "link", se.fit = TRUE)
gam_ci_df <- data.frame(
  fit = gam_link_ci$fit,
  lwr = gam_link_ci$fit - 1.96 * gam_link_ci$se.fit,
  upr = gam_link_ci$fit + 1.96 * gam_link_ci$se.fit
)

# Apply plogis and rename columns before adding to model_outputs_2
logit_lm_ci_df <- as.data.frame(logit_lm_ci) %>%
  mutate(across(everything(), plogis)) %>%
  rename(
    logit_lm_fit = fit,
    logit_lm_lwr = lwr,
    logit_lm_upr = upr
  )

# Do the same for the other models
base_lm_ci_df <- as.data.frame(base_lm_ci) %>%
  mutate(across(everything(), plogis)) %>%
  rename(
    base_lm_fit = fit,
    base_lm_lwr = lwr,
    base_lm_upr = upr
  )

glm_ci_df <- glm_ci_df %>%
  mutate(across(everything(), plogis)) %>%
  rename(
    glm_fit = fit,
    glm_lwr = lwr,
    glm_upr = upr
  )

gam_ci_df <- gam_ci_df %>%
  mutate(across(everything(), plogis)) %>%
  rename(
    gam_fit = fit,
    gam_lwr = lwr,
    gam_upr = upr
  )

# Bind all columns together with KK_1_EPG into model_outputs_2
model_outputs_2 <- bind_cols(
  combined_df_2["KK_1_EPG"],
  logit_lm_ci_df,
  base_lm_ci_df,
  glm_ci_df,
  gam_ci_df
)

# ---- Combine into summary table ----
summary_table_2 <- bind_rows(
  get_prediction_summary(model_outputs_2, "logit_lm_fit", "logit_lm_lwr", "logit_lm_upr"),
  get_prediction_summary(model_outputs_2, "base_lm_fit", "base_lm_lwr", "base_lm_upr"),
  get_prediction_summary(model_outputs_2, "glm_fit", "glm_lwr", "glm_upr"),
  get_prediction_summary(model_outputs_2, "gam_fit", "gam_lwr", "gam_upr")
) %>%
  mutate(Model = c("Logit_lm", "Base_lm", "GLM", "GAM")) %>%
  select(Model, everything())

# ---- View result ----
print(summary_table_2)

####Combined_df_3

# ---- Create a new data frame with only the KK_1_EPG column ----
model_outputs_3 <- combined_df_3 %>%
  select(KK_1_EPG)

# ---- Predict with confidence intervals for each model ----
logit_lm_ci <- predict(model_logit_lm, newdata = model_outputs_3, interval = "confidence")
base_lm_ci  <- predict(model_logit, newdata = model_outputs_3, interval = "confidence")

# GLM
glm_link_ci <- predict(model_glm, newdata = model_outputs_3, type = "link", se.fit = TRUE)
glm_ci_df <- data.frame(
  fit = glm_link_ci$fit,
  lwr = glm_link_ci$fit - 1.96 * glm_link_ci$se.fit,
  upr = glm_link_ci$fit + 1.96 * glm_link_ci$se.fit
)

# GAM
gam_link_ci <- predict(model_gam, newdata = model_outputs_3, type = "link", se.fit = TRUE)
gam_ci_df <- data.frame(
  fit = gam_link_ci$fit,
  lwr = gam_link_ci$fit - 1.96 * gam_link_ci$se.fit,
  upr = gam_link_ci$fit + 1.96 * gam_link_ci$se.fit
)

# Apply plogis and rename columns
logit_lm_ci_df <- as.data.frame(logit_lm_ci) %>%
  mutate(across(everything(), plogis)) %>%
  rename(
    logit_lm_fit = fit,
    logit_lm_lwr = lwr,
    logit_lm_upr = upr
  )

base_lm_ci_df <- as.data.frame(base_lm_ci) %>%
  mutate(across(everything(), plogis)) %>%
  rename(
    base_lm_fit = fit,
    base_lm_lwr = lwr,
    base_lm_upr = upr
  )

glm_ci_df <- glm_ci_df %>%
  mutate(across(everything(), plogis)) %>%
  rename(
    glm_fit = fit,
    glm_lwr = lwr,
    glm_upr = upr
  )

gam_ci_df <- gam_ci_df %>%
  mutate(across(everything(), plogis)) %>%
  rename(
    gam_fit = fit,
    gam_lwr = lwr,
    gam_upr = upr
  )

# Combine with KK_1_EPG into model_outputs_3
model_outputs_3 <- bind_cols(
  combined_df_3["KK_1_EPG"],
  logit_lm_ci_df,
  base_lm_ci_df,
  glm_ci_df,
  gam_ci_df
)

# ---- Combine into summary table ----
summary_table_3 <- bind_rows(
  get_prediction_summary(model_outputs_3, "logit_lm_fit", "logit_lm_lwr", "logit_lm_upr"),
  get_prediction_summary(model_outputs_3, "base_lm_fit", "base_lm_lwr", "base_lm_upr"),
  get_prediction_summary(model_outputs_3, "glm_fit", "glm_lwr", "glm_upr"),
  get_prediction_summary(model_outputs_3, "gam_fit", "gam_lwr", "gam_upr")
) %>%
  mutate(Model = c("Logit_lm", "Base_lm", "GLM", "GAM")) %>%
  select(Model, everything())

# ---- View result ----
print(summary_table_3)
print(summary_table_2)
