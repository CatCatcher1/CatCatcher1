if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiomeMarker")
remotes::install_github("biobakery/Maaslin2")

library(phyloseq)
library(Maaslin2)
library(microbiomeMarker)

otu_table = otu_table(ps)
sample_data = sample_data(ps)
demo_output = "C:/Users/Batuhan/Desktop/Medipol Staj/new_data/MaAsLin2"

library(Maaslin2)
library(Maaslin2)
input_data <- system.file(
  'extdata',otu_table_ , package="Maaslin2")
input_metadata <-system.file(
  'extdata',sample_data_, package="Maaslin2")
fit_data <- Maaslin2(
  input_data, input_metadata, "demo_output",
  fixed_effects = c('Genus',"Family"),
  random_effects = c('site', 'subject'),
  standardize = FALSE)

###
library(Maaslin2)
otu_table <- otu_table(ps)
sample_data <- sample_data(ps)
demo_output <- "C:/Users/Batuhan/Desktop/Medipol Staj/new_data/MaAsLin2"

fit_data <- Maaslin2(
  input_data = otu_table,
  input_metadata = sample_data,
  output = demo_output,
  fixed_effects = NULL, 
  random_effects = NULL, 
  standardize = FALSE
)

###

library(Maaslin2)
input_data <- system.file(
  'extdata',as.data.frame(otu_table(ps)), package="Maaslin2")
input_metadata <-system.file(
  'extdata',as.data.frame(sample_data(ps)), package="Maaslin2")
fit_data <- Maaslin2(
  input_data, input_metadata, 'C:/Users/Batuhan/Desktop/Medipol Staj/new_data/MaAsLin2',
  fixed_effects = c('Subject', 'Gender', 'Day'),
  random_effects = c('When'),
  standardize = FALSE)

###

library(Maaslin2)

# Assuming ps is a matrix containing your OTU data
otu_table <- otu_table(ps)  # Ensure otu_table is a matrix

sample_data <- sample_data(ps)  # Check the output type here

otu_table_frame <- as.data.frame(otu_table)
sample_data_frame <- as.data.frame(sample_data)

demo_output <- "C:/Users/Batuhan/Desktop/Medipol Staj/new_data/MaAsLin2"

fit_data <- Maaslin2(
  input_data = otu_table_frame,
  input_metadata = sample_data_frame,
  output = demo_output,
  fixed_effects = c('Sample', 'BodySite'),  # No fixed effects
  random_effects = NULL, 
  standardize = FALSE
)
###

fit_data2 = Maaslin2(
  input_data = otu_table, 
  input_metadata = sample_data, 
  output = "demo_output2", 
  fixed_effects = c("Sample", "BodySite"))

###
install.packages("devtools")
library("devtools")
install_github("biobakery/Maaslin2")
library(Maaslin2)
fit_data_taxa_output <- Maaslin2(input_data = data.frame(otu_table(ps.prop.genus)), 
                                 input_metadata = data.frame(sample_data(ps.prop.genus)),
                                 output = "C:/Users/Batuhan/Desktop/Medipol Staj/new_data/MaAsLin2_v2", 
                                 fixed_effects = c("BodySite"), min_prevalence = 0.10, 
                                 normalization = "TSS", transform="LOG", analysis_method = "LM", standardize = FALSE, 
                                 plot_heatmap = TRUE)
