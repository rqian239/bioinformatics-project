library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
library(scales)
library(reshape2)

expression_df <- readr::read_tsv("./results/SRP073813-HUGO-cleaned.tsv") %>% 
  mutate(first_mapped_hugo = make.unique(as.character(first_mapped_hugo))) %>%
  tibble::column_to_rownames("first_mapped_hugo")

expression_df$med_expr <- apply(expression_df, 1, FUN = median, na.rm = T)

ggplot(expression_df, aes(x = expression_df$med_expr)) +
  geom_density(fill = "#69b3a2", color = "#e9ecef", alpha = 0.8) +
  xlab("log_2(count)") + 
  ggtitle("Density Plot for Median Expression Counts")
  
ggsave("med_expr_density.png")


ggplot(expression_df, aes(x = expression_df$med_expr)) +
  geom_density(fill = "#69b3a2", color = "#e9ecef", alpha = 0.8) +
  xlim(0, 20) +
  xlab("log_2(count)") + 
  ggtitle("Density Plot for Median Expression Counts (Bounded Before 20)")

ggsave("med_expr_density_bounded.png")

