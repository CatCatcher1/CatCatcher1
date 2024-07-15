install.packages("csv")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages("tidyverse")
library(tidyverse)
library(phyloseq)
library(vegan)
library(csv)
library(ape)
library(devtools)
library(pairwiseAdonis)
###
set.seed(81)
###
head(sample_data(ps))
###
bray <- ordinate(
  physeq = ps,
  method = "PCoA", 
  distance = "bray" 
)
###
plot_ordination(
  physeq = ps,                                                         
  ordination = bray)+                                               
  geom_point(aes(fill = Sample, shape = BodySite), size = 6) 
###
colors <- c("lightblue", "maroon", "white", "aquamarine3", "bisque3")

plot_ordination(
  physeq = ps,                                                          
  ordination = bray)+                                                       
  geom_point(aes(fill = BodySite, shape = BodySite), size = 4) +                         
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(21, 22))+
  theme_linedraw() +                                                         
  theme(                             
    legend.title = element_blank(),                                         
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))
###
random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
plot(random_tree)
###
ps2 = merge_phyloseq(ps, random_tree)
ps2
###
uni <- ordinate(
  physeq = ps2, 
  method = "PCoA", 
  distance = "unifrac"
)

colors <- c("lightblue", "maroon")

plot_ordination(
  physeq = ps2,                                                          
  ordination = uni)+                                                
  geom_point(aes(fill = BodySite, shape = BodySite), size = 6) +                         
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  theme_linedraw() +                                                      
  theme(                             
    legend.title = element_blank(),                                      
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21))) 