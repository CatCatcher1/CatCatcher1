library(phyloseq)
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)

richness <- estimate_richness(ps)
head(richness)
###
plot_richness(ps)
###
plot_richness(ps, x="BodySite", measures="Shannon", color = "BodySite")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

hist(richness$Shannon, main="Shannon index", xlab="")
###
a_my_comparisons <- list( c("Fecal", "Saliva"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c( "*", "ns"))

plot_richness(ps, x="BodySite", measures="Shannon", color = "BodySite")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)