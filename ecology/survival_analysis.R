# Load necessary libraries
library(survminer)
require(survival)

# Load the lung cancer dataset
?lung # Read about lung cancer dataset
data(lung)

# Perform Kaplan-Meier survival analysis and visualization by sex
df <- survfit(Surv(time, status) ~ sex, data = lung)
ggsurvplot(
  df,
  ggtheme = theme_minimal(),
  linetype = "strata",
  risk.table = TRUE,
  risk.table.col = "strata",
  conf.int = TRUE,
  pval = TRUE,
  palette = c("#fdbf11", "#12719e"),
  legend.labs = c("Male", "Female")
)

# Perform Cox proportional hazards regression modeling
res.cox <- coxph(Surv(time, status) ~ ph.karno * age, data = lung)

# Display summary of Cox regression results
summary(res.cox, conf.int = FALSE)

# Generate forest plot for Cox regression
ggforest(res.cox, data = lung)

# Create interaction term ph.karno_age
lung$ph.karno_age <- lung$ph.karno * lung$age

# Perform Cox proportional hazards regression with interaction term
res.cox2 <- coxph(Surv(time, status) ~ ph.karno + age + ph.karno_age, data = lung)

# Display summary of second Cox regression model
summary(res.cox2 , conf.int = FALSE)

# Generate forest plot for second Cox regression model
ggforest(res.cox2, data = lung)
