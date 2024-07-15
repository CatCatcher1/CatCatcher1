install.packages("devtools")
library("devtools")
install_github("biobakery/Maaslin2")
library(Maaslin2)
fit_data_taxa_output <- Maaslin2(input_data = data.frame(otu_table(ps.prop.genus)), 
                                 input_metadata = data.frame(sample_data(ps.prop.genus)),
                                 output = "C:/Users/Batuhan/Desktop/Medipol Staj/new_data/MaAsLin2_v2", 
                                 fixed_effects = "BodySite", min_prevalence = 0.10, 
                                 normalization = "TSS", transform="LOG", analysis_method = "LM", standardize = FALSE, 
                                 plot_heatmap = TRUE)