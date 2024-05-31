# Learn to reshape data into tidy data format for tidyverse and vice versa

# Example data
ID <- c(1, 2, 3)
Name <- c("Andika", "Angga", "Joko")
Variable1 <- c(10, 20, 30)
Variable2 <- c(15, 25, 35)
wide_data <- data.frame(ID, Name, Variable1, Variable2)

# Using reshape2 package
library(reshape2)
# Melting wide-format data into long-format
melted_data <- melt(wide_data, id.vars = c("ID", "Name"), measure.vars = c("Variable1", "Variable2"))
# Casting long-format data into wide-format
reshaped_data <- dcast(melted_data, ID + Name ~ variable)

# Using tidyr package
library(tidyr)
# Pivoting longer data into wider format
tidyr_wide_data <- pivot_wider(melted_data, names_from = "variable", values_from = "value")
# Pivoting wider data into longer format
tidyr_long_data <- pivot_longer(wide_data, cols = -c(ID, Name), names_to = "Variable", values_to = "Value")

# Compare the results

print(reshaped_data)

print(tidyr_wide_data)

print(tidyr_long_data)
