 # install.packages("nlme")
library(nlme)
library(lme4)
library(readxl)
library(ggplot2)



ARM <- read_excel("ARM.xlsx")
linearmixed = ARM

#preparing data
linearmixed = as.data.frame(linearmixed)
rownames(linearmixed) = linearmixed[,1]
linearmixed = linearmixed[,-1]


# Define the model formula
model_formula <- "pseudotime ~ amyloid_pathology + tau_pathology + APOE_genotype + TREM2_genotype + (1 | sample_ID)"

# Fit the linear mixed-effects model
model <- lmer(model_formula, data = linearmixed)

# Print the model summary
summary(model)
# Extract fixed effects from the model
fixed_effects <- fixef(model)

fixed_effects_df <- data.frame(Term = names(fixed_effects), Estimate = fixed_effects)

# Scatter plot of estimates vs terms
ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 3)), vjust = -0.5) +  # Adjust vjust for vertical positioning
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a line from 0
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0))  # Adjust the hjust value for x-axis labels



# Filter out the intercept term
fixed_effects_df <- fixed_effects_df[fixed_effects_df$Term != "(Intercept)", ]

ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 3)), hjust = -0.2) +  # Add labels for estimates
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a line from 0
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))

# other way to visaulize

fixed_effects_df <- fixed_effects_df[fixed_effects_df$Term != "(Intercept)", ]

# Create the plot
ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 4)), hjust = -0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1)) +
  scale_x_continuous(
    limits = c(-0.1, 0.1),
    breaks = seq(-0.1, 0.1, by = 0.04),  
    labels = scales::number_format(accuracy = 0.01)
  )


#for Motile group

Motile <- read_excel("Motile.xlsx")
linearmixed = Motile

linearmixed = as.data.frame(linearmixed)
rownames(linearmixed) = linearmixed[,1]
linearmixed = linearmixed[,-1]



# Define the model formula
model_formula <- "pseudotime ~ amyloid_pathology + tau_pathology + APOE_genotype + TREM2_genotype + (1 | sample_ID)"

# Fit the linear mixed-effects model
model <- lmer(model_formula, data = linearmixed)

# Print the model summary
summary(model)
# Extract fixed effects from the model
fixed_effects <- fixef(model)

fixed_effects_df <- data.frame(Term = names(fixed_effects), Estimate = fixed_effects)

# Scatter plot of estimates vs terms
library(ggplot2)

ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 3)), vjust = -0.5) +  # Adjust vjust for vertical positioning
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a line from 0
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0))  # Adjust the hjust value for x-axis labels


library(ggplot2)

# Filter out the intercept term
fixed_effects_df <- fixed_effects_df[fixed_effects_df$Term != "(Intercept)", ]

ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 3)), hjust = -0.2) +  # Add labels for estimates
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a line from 0
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))

# other way again with adjusting some axis points
fixed_effects_df <- fixed_effects_df[fixed_effects_df$Term != "(Intercept)", ]

# Create the plot
ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 4)), hjust = -0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1)) +
  scale_x_continuous(
    limits = c(-0.1, 0.2),
    breaks = seq(-0.1, 0.2, by = 0.04),  
    labels = scales::number_format(accuracy = 0.01)
  )


# again for ARM group but here I also consider APOE isoforms seperately

ARM_with_APOEgenotype <- read_excel("ARM_with_APOEgenotype.xlsx")
linearmixed2 = ARM_with_APOEgenotype

linearmixed2 = as.data.frame(linearmixed2)
rownames(linearmixed2) = linearmixed2[,1]
linearmixed2 = linearmixed2[,-1]


# Define the model formula with APOE isoforms
model_formula <- model_formula <- "pseudotime ~ amyloid_pathology + tau_pathology + `APOE_genotype[E4/E4]` + `APOE_genotype[E3/E4]` + TREM2_genotype + (1 | sample_ID)"

# Fit the linear mixed-effects model
model <- lmer(model_formula, data = linearmixed2)

# Print the model summary
summary(model)
# Extract fixed effects from the model
fixed_effects <- fixef(model)

fixed_effects_df <- data.frame(Term = names(fixed_effects), Estimate = fixed_effects)

# Scatter plot of estimates vs terms

ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 3)), vjust = -0.5) +  # Adjust vjust for vertical positioning
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a line from 0
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0))  # Adjust the hjust value for x-axis labels



# Filter out the intercept term
fixed_effects_df <- fixed_effects_df[fixed_effects_df$Term != "(Intercept)", ]

ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 4)), hjust = -0.2) +  # Add labels for estimates
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a line from 0
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))



fixed_effects_df <- fixed_effects_df[fixed_effects_df$Term != "(Intercept)", ]

# Create the plot
ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 4)), hjust = -0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1)) +
  scale_x_continuous(
    limits = c(-0.1, 0.2),
    breaks = seq(-0.1, 0.2, by = 0.04),  
    labels = scales::number_format(accuracy = 0.01)
  )



# again for Motile group but here I also consider APOE isoforms seperately

Motile_with_APOEgenotype_new <- read_excel("Motile_with_APOEgenotype.xlsx")

linearmixed2 = Motile_with_APOEgenotype_new

linearmixed2 = as.data.frame(linearmixed2)
rownames(linearmixed2) = linearmixed2[,1]
linearmixed2 = linearmixed2[,-1]


# Define the model formula
model_formula <- model_formula <- "pseudotime ~ amyloid_pathology + tau_pathology + `APOE_genotype[E4/E4]` + `APOE_genotype[E3/E4]` + TREM2_genotype + (1 | sample_ID)"

# Fit the linear mixed-effects model
model <- lmer(model_formula, data = linearmixed2)

# Print the model summary
summary(model)
# Extract fixed effects from the model
fixed_effects <- fixef(model)

fixed_effects_df <- data.frame(Term = names(fixed_effects), Estimate = fixed_effects)

# Scatter plot of estimates vs terms

ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 3)), vjust = -0.5) +  # Adjust vjust for vertical positioning
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a line from 0
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0))  # Adjust the hjust value for x-axis labels


##
# Filter out the intercept term
fixed_effects_df <- fixed_effects_df[fixed_effects_df$Term != "(Intercept)", ]

ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 4)), hjust = -0.2) +  # Add labels for estimates
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a line from 0
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))


###

fixed_effects_df <- fixed_effects_df[fixed_effects_df$Term != "(Intercept)", ]

# Create the plot
ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(position = position_dodge(0.5), size = 2, color = "skyblue") +
  geom_text(aes(label = round(Estimate, 4)), hjust = -0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Linear Mixed Effects Model Coefficients",
       x = "Estimate", y = "Fixed Effects") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1)) +
  scale_x_continuous(
    limits = c(-0.1, 0.2),
    breaks = seq(-0.1, 0.2, by = 0.04),  
    labels = scales::number_format(accuracy = 0.01)
  )



