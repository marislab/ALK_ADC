# script to generate the manuscript figures 1A-B
# Ensure the working directory is set to the directory with the code and data available

library(tidyverse)

# 1. Normal Tissue Plot -------------
ot_data <- readRDS('ot_data_expr.rds')
histologies.subset.nb <- read.delim('histologies.txt', header = T)

# generate plot
gtex.plot <- ot_data %>% 
  rownames_to_column(var = 'gene') %>% 
  gather(key = 'Kids_First_Biospecimen_ID', value = 'TPM', -gene) %>% 
  inner_join(., histologies.subset.nb, by = 'Kids_First_Biospecimen_ID') %>% 
  #mutate(TPM = TPM + 0.01) %>%
  mutate(log2TPM = log2(TPM + 1)) %>% 
  mutate(cohort = factor(cohort, levels = c('TARGET', 'GMKF', 'GTEx'))) %>%
  ggplot(., aes(reorder(histology, log2TPM, median), log2TPM)) +
  stat_boxplot(geom ='errorbar', width = 0.2, linewidth = 1) +
  geom_boxplot(lwd = 1, fatten = 0.7, width = 0.5, outlier.size = 2, outlier.alpha = 0, fill = 'lightgray') + # used outlier.alpha = 0.3
  theme_classic() + 
  #scale_fill_manual(values = c('darkgoldenrod2', 'deepskyblue3', 'lightgrey')) +
  coord_flip() +
  labs(title = 'ALK expression in normal tissues', x = '', y = 'log2(TPM + 1)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "right",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 25),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        text=element_text(face="bold"))


# 2. ALK Outlier Expression Plot ------------------
treehouse_modified <- readRDS('2025-06-20_manuscript_figures_code/ALK_expression_outlier_histologies.rds')

createPlots <- function(disease_vector){
  outliers <- treehouse_modified[treehouse_modified$disease %in% disease_vector,]
  
  # make treehouse_modified histologies in sentence case
  outliers$disease <- str_to_title(outliers$disease)
  
  # plot
  outlier_plot <- outliers |>
    mutate(log2TPM = log2(TPM + 1)) |> 
    ggplot(aes(x=reorder(disease, log2TPM, FUN=median), y=log2TPM)) + 
    stat_boxplot(geom ='errorbar', width = 0.2, linewidth = 1) +
    geom_boxplot(lwd = 1, fatten = 0.7, width = 0.5, outlier.size = 5, outlier.alpha = 0, fill = 'lightgray') +
    scale_y_continuous(breaks = round(seq(0, max(outliers$TPM), by = 5))) +
    theme_classic() + 
    coord_flip() +
    labs(title = 'Outliers in ALK expression', x = '', y = 'log2(TPM + 1)') + 
    theme(plot.title = element_text(hjust = 0.5, size = 25),
      legend.position = "none",
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 25),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
      text=element_text(face="bold"))
  
  
  remaining <- treehouse_modified
  # make treehouse_modified histologies in sentence case
  remaining$disease <- str_to_title(remaining$disease)
  
  remaining_plot <- remaining |>
    mutate(log2TPM = log2(TPM + 1)) |> 
    ggplot(aes(x=reorder(disease, log2TPM, FUN = median), y=log2TPM)) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, width = 0.5, outlier.size = 3, outlier.alpha = 0, fill = 'lightgray') +
    scale_y_continuous(breaks = round(seq(0, max(outliers$TPM), by = 5))) +
    theme_classic() + 
    coord_flip() +
    labs(title = 'ALK expression in all histologies', x = '') + 
    theme(plot.title = element_text(hjust = 0.5, size = 15),
      legend.position = "none",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))
  
  final <- list(outlier_plot, remaining_plot)
  return(final)
  
}


# median TPM in 95 percentile & %tage of Samples in 99th percentile
# another strategy would be to find out %tage of samples with expression in 90th or 95th percentile
# which diseases are in and above 95 percentile of median TPM
median_TPMs <- treehouse_modified %>% 
  group_by(disease) %>% 
  mutate(median_TPM = median(TPM)) %>% 
  dplyr::select(disease, median_TPM) %>% 
  distinct()

quantile(median_TPMs$median_TPM, c(.90, .95, .99))
# 90%       95%       99% 
# 2.760098  6.555000 18.084952 

top95percentile <- median_TPMs %>% 
  filter(median_TPM >= 6.555000) %>% 
  .$disease


# %tage of samples in each disease to be above 99 percentile
quantile(treehouse_modified$TPM, c(.90, .95, .99)) 
# 90%       95%       99% 
# 2.790143  5.899500 31.040332 

counts_99percentile <- treehouse_modified %>% 
  group_by(disease) %>% 
  mutate(count_samples_99 = sum(TPM >= 31.040332)) %>% 
  dplyr::select(disease, count_samples_99) %>% 
  distinct() %>% 
  filter(count_samples_99 > 1)



# take a union of two sets
median95_counts99 <- unique(union(top95percentile, counts_99percentile$disease))
median95_counts99 <- median95_counts99[-5]
final <- createPlots(disease_vector = median95_counts99)


p1 <- final[[1]]
p2 <- gtex.plot
final_plot <- gridExtra::grid.arrange(p1, p2, ncol = 2)

ggsave(final_plot, filename = paste0(Sys.Date(),'_Figure1_A_B_plot.pdf'), width = 20, height = 12)
