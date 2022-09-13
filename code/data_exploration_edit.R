# Data exploration
  ## look at distribution and abundance of S. musiva among trees (violin plots across disease scores)
  ## visualizing differences between healthy, non-canker and canker samples (ordination - from prep script)
  ## visualizing differences between genotypes (east/west comparison)
  ## correlation of S. musiva abundance/presence with other OTUs - vectored NMDS plots, pairwise scatterplots
  ## visualizing ecological groups across trees (get funguild, then stacked bar plots)
# Using 97% clustered only
# Using samples with minDepth of 500 reads

library(tidyverse)
library(phyloseq)
library(vegan)
library(ggthemes)
library(dada2)
library(ggplot2)
library(ggpubr)
library(indicspecies)
library(viridis)
library(pairwiseAdonis)
library(magrittr)
library(reshape2)
#library(devtools)

sessionInfo()

#### Preparing the data ####
depth <- read.csv("../output/analysis/May_data(outdated)/depth.csv")
clustered.meta <- read.csv("../output/analysis/June_data/clustered.meta.csv", row.names = 1)
clustered.phy <- readRDS("../output/analysis/June_data/clustered.phy.rds")
sample_data(clustered.phy) <- clustered.meta
clust.rich.phy <- readRDS("../output/analysis/April_data_(outdated)/clust.rich.phy.rds") #this phyloseq object has proportional read counts, and all OTUs are included 
clust.prop.phy <- clustered.phy %>% transform_sample_counts(function(x){x*min(sample_sums(.)/sum(x))})

minDepth <- 500
clust2.5.phy <- (clust.prop.phy %>% filter_taxa(function(x) { sum(x>0) > 0.025*ntaxa(.) }, TRUE) %>%
                   prune_samples(sample_sums(.)>minDepth,.) %>% 
                   filter_taxa(function(x) {sum(x) > 0}, TRUE))

clust.otu <- as.data.frame(otu_table(clust.prop.phy))

healthy.phy <- subset_samples(clust.prop.phy, Disease=="Healthy")
healthy.2.5.phy <- (healthy.phy %>% filter_taxa(function(x) { sum(x>0) > 0.025*ntaxa(.) }, TRUE) %>%
                   prune_samples(sample_sums(.)>minDepth,.) %>% 
                   filter_taxa(function(x) {sum(x) > 0}, TRUE))
clust.dis.phy <- subset_samples(clust.prop.phy, Disease=="Canker"|Disease=="Non-canker")
clust.dis.2.5.phy <- (healthy.phy %>% filter_taxa(function(x) { sum(x>0) > 0.025*ntaxa(.) }, TRUE) %>%
                     prune_samples(sample_sums(.)>minDepth,.) %>% 
                     filter_taxa(function(x) {sum(x) > 0}, TRUE))
canker.phy <- subset_samples(clust.prop.phy, Disease=="Canker")
canker.2.5.phy <- (canker.phy %>% filter_taxa(function(x) { sum(x>0) > 0.025*ntaxa(.) }, TRUE) %>%
                        prune_samples(sample_sums(.)>minDepth,.) %>% 
                        filter_taxa(function(x) {sum(x) > 0}, TRUE))
noncanker.phy <- subset_samples(clust.prop.phy, Disease=="Non-canker")
noncanker.2.5.phy <- (noncanker.phy %>% filter_taxa(function(x) { sum(x>0) > 0.025*ntaxa(.) }, TRUE) %>%
                     prune_samples(sample_sums(.)>minDepth,.) %>% 
                     filter_taxa(function(x) {sum(x) > 0}, TRUE))
pop.phy <- subset_samples(clust.prop.phy, Pop=="Core"|Pop=="South"|Pop=="BC"|Pop=="Columbia")
pop2.5.phy <- subset_samples(clust2.5.phy, Pop=="Core"|Pop=="South"|Pop=="BC"|Pop=="Columbia")
#pop.20.phy <- subset_samples(clust.20.phy, Pop=="Core"|Pop=="South"|Pop=="BC"|Pop=="Columbia")
core.s.phy <- subset_samples(clust.prop.phy, Pop=="Core"|Pop=="South")

saveRDS(clust.prop.phy, "../output/analysis/June_data/clust.prop.phy.rds")
otu_table(clust.prop.phy) %>% write.csv("../output/analysis/June_data/OTU.prop.table.csv")
tax_table(clust.prop.phy) %>% write.csv("../output/analysis/June_data/OTU.prop.taxonomy.table.csv")

saveRDS(healthy.2.5.phy, "../output/analysis/June_data/healthy.2.5.phy.rds")
otu_table(healthy.2.5.phy) %>% write.csv("../output/analysis/June_data/OTU.prop.table.csv")
tax_table(healthy.2.5.phy) %>% write.csv("../output/analysis/June_data/OTU.health2.5.taxonomy.table.csv")

saveRDS(noncanker.2.5.phy, "../output/analysis/June_data/noncanker.2.5.phy.rds")
otu_table(noncanker.2.5.phy) %>% write.csv("../output/analysis/June_data/OTU.prop.table.csv")
tax_table(noncanker.2.5.phy) %>% write.csv("../output/analysis/June_data/OTU.noncank2.5.taxonomy.table.csv")

saveRDS(canker.2.5.phy, "../output/analysis/June_data/canker.2.5.phy.rds")
otu_table(canker.2.5.phy) %>% write.csv("../output/analysis/June_data/OTU.prop.table.csv")
tax_table(canker.2.5.phy) %>% write.csv("../output/analysis/June_data/OTU.cank2.5.taxonomy.table.csv")


#saveRDS(clust2.5.phy, "output/analysis/June_data/clust2.5.phy.rds")
#otu_table(clust2.5.phy) %>% write.csv("output/analysis/June_data/OTU2.5.table.csv")
#tax_table(clust2.5.phy) %>% write.csv("output/analysis/June_data/OTU2.5.taxonomy.table.csv")

#saveRDS(clust.20.phy, "output/analysis/June_data/clust.20.phy.rds")
#otu_table(clust.20.phy) %>% write.csv("output/analysis/June_data/OTU.20.table.csv")
#tax_table(clust.20.phy) %>% write.csv("output/analysis/June_data/OTU.20.taxonomy.table.csv")

#### Checking data quality ####
depth.hist <- ggplot(depth, aes(x=Depth)) +
  geom_histogram(binwidth = 1000)
accum.curve <- specaccum(clust.otu)
plot(accum.curve, add = FALSE, random = FALSE, ci =2, ci.type = "bar", ci.col = 'red', lty = 1)
OTU.depth.curve <- ggplot(depth, aes(x=Depth, y=OTU_count)) +
  geom_point() +
  geom_smooth(method = "loess")
# Posy had an issue with using this formula
plot(OTU.depth.curve, xlim=c(0,6000), ylim=c(0,40))

#### Taxonomic distribution ####
theme_set(theme_bw())
clust.prop.phy@tax_table@.Data %<>% parse_taxonomy_default()
pop.phy@tax_table@.Data %<>% parse_taxonomy_greengenes()
core.s.phy@tax_table@.Data %<>% parse_taxonomy_greengenes()

orders <- plot_bar(clust2.5.phy, x="Sample", fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack", color="gray33", size = 0.05) +
  facet_wrap(~Disease, nrow = 1, scales = "free_x", drop = TRUE)+ 
  theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), axis.text.x = element_blank()) + 
  ylab("Proportional abundance") +
  scale_fill_viridis(begin = 0, end = 1, discrete = T, option = "C") 

healthy.2.5.order <- plot_bar(healthy.2.5.phy, x="Sample", fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack", color="gray33", size = 0.05) +
  facet_wrap(~Disease, nrow = 1, scales = "free_x", drop = TRUE)+ 
  theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), axis.text.x = element_blank()) + 
  ylab("Proportional abundance") +
  scale_fill_viridis(begin = 0, end = 1, discrete = T, option = "C")

healthy.order <- plot_bar(healthy.phy, x="Sample", fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack", color="gray33", size = 0.05) +
  facet_wrap(~Disease, nrow = 1, scales = "free_x", drop = TRUE)+ 
  theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), axis.text.x = element_blank()) + 
  ylab("Proportional abundance") +
  scale_fill_viridis(begin = 0, end = 1, discrete = T, option = "C")

pop.orders <- plot_bar(pop.phy, x="Sample", fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack", color="gray33") +
  facet_wrap(~Pop, nrow = 1, scales = "free_x", drop = TRUE)+ 
  theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), axis.text.x = element_blank()) + 
  ylab("Proportional abundance") +
  scale_fill_viridis(begin = 0, end = 1, discrete = T, option = "C")

core.s.orders <- plot_bar(core.s.phy, x="Sample", fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack", color="gray33") +
  facet_wrap(~Pop, nrow = 1, scales = "free_x", drop = TRUE)+ 
  theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), axis.text.x = element_blank()) + 
  ylab("Proportional abundance") +
  scale_fill_viridis(begin = 0, end = 1, discrete = T, option = "C")


# *****************  core.tax.phy not defined

# core.tax.phy <- subset_taxa(clust.prop.phy, fill=='core')
# core.tax.phy <- prune_samples(sample_sums(core.tax.phy)>0, core.tax.phy)
# core.tax.phy@tax_table@.Data %<>% parse_taxonomy_greengenes()
# core.tax.gen <- plot_bar(core.tax.phy, x="Sample", fill = "Genus") + 
#   geom_bar(aes(fill=Genus), stat="identity", position="stack", color="gray33") +
#   facet_wrap(~Disease, nrow = 1, scales = "free_x", drop = TRUE)+ 
#   theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), axis.text.x = element_blank()) + 
#   ylab("Proportional abundance") +
#   scale_fill_viridis(begin = 0, end = 1, discrete = T, option = "C")

# *****************   clust.rich.phy not defined

#### Calculating and testing alpha diversity ####
otu.rich <- otu_table(clust.rich.phy)
meta.rich <- sample_data(clust.rich.phy)
otu.shannon <- diversity(otu.rich, index = "shannon", MARGIN = 1, base = exp(1))
otu.shannon.df <- as.data.frame(otu.shannon)
otu.shannon.df$Disease <- meta.rich$Disease
# shannon diversity boxplot by population
otu.shannon.box <- ggplot(otu.shannon.df, aes(x=Disease, y=otu.shannon, fill = Disease)) +
  geom_boxplot() +
  ylab("Shannon diversity") +
  theme(legend.position = "none")
# ANOVA of alpha diversity by disease group
shannon.aov <- aov(otu.shannon ~ Disease, data = otu.shannon.df)
summary(shannon.aov)
TukeyHSD(shannon.aov)

# merge all 2.5 phyloseq objects and perform same diversity tests
merged.2.5.phy <- merge_phyloseq(healthy.2.5.phy,noncanker.2.5.phy,canker.2.5.phy)
merged.rich <- otu_table(merged.2.5.phy)
meta.merged <- sample_data(merged.2.5.phy)
merged.shannon <- diversity(merged.rich, index = "shannon", MARGIN = 1, base = exp(1))
merged.shannon.df <- as.data.frame(merged.shannon)
merged.shannon.df$Disease <- meta.merged$Disease
merged.shannon.box <- ggplot(merged.shannon.df, aes(x=Disease, y=merged.shannon, fill = Disease)) +
    geom_boxplot() +
    ylab("Shannon diversity") +
    theme(legend.position = "none")
shannon.merged.aov <- aov(merged.shannon ~ Disease, data = merged.shannon.df)
summary(shannon.merged.aov)
TukeyHSD(shannon.merged.aov)

#### Ordinations ####
# all trees at 2.5% prevalence
clust.bray <- clust2.5.phy %>% phyloseq::distance("bray") %>% sqrt
clust.nmds <- metaMDS(clust.bray, trymax = 200, parallel = 10)
  ## no convergence, run 200 stress 0.198
  ## outliers (below) removed: no convergence, run 200 stress=0.2098
clust.nmds.dat <- scores(clust.nmds) %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(clust2.5.phy) %>% data.frame %>% rownames_to_column("sampID"))
ggplot(clust.nmds.dat,aes(x=NMDS1,y=NMDS2,fill=Disease,label = sample_names(clust2.5.phy))) +
  geom_point(shape=21,size=3)+
  geom_text() +
  stat_ellipse() +
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few()


merged.bray <- merged.2.5.phy %>% phyloseq::distance("bray") %>% sqrt
merged.nmds <- metaMDS(merged.bray, trymax = 200, parallel = 10)
## no convergence, run 200 stress 0.198
## outliers (below) removed: no convergence, run 200 stress=0.2098
merged.nmds.dat <- scores(merged.nmds) %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(merged.2.5.phy) %>% data.frame %>% rownames_to_column("sampID"))
ggplot(merged.nmds.dat,aes(x=NMDS1,y=NMDS2,fill=Disease,label = sample_names(merged.2.5.phy))) +
  geom_point(shape=21,size=3)+
  geom_text() +
  stat_ellipse() +
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few()


  ## outliers: myco.03.02.f, myco.06.04.d, myco.03.04.a, myco.07.02.c, myco.07.05.d, myco.04.02.e
# removing outliers
clust.in.phy <- subset_samples(clust2.5.phy, sample_names(clust2.5.phy) != c("myco.03.02.f", "myco.06.04.d", 
"myco.03.04.a")) 
clust.in.phy <- subset_samples(clust.in.phy, sample_names(clust.in.phy) != c("myco.07.02.c")) 
clust.in.phy <- subset_samples(clust.in.phy, sample_names(clust.in.phy) != c("myco.07.05.d"))
clust.in.phy <- subset_samples(clust.in.phy, sample_names(clust.in.phy) != c("myco.04.02.e"))
clust.in.phy <- subset_samples(clust.in.phy, sample_names(clust.in.phy) != c("myco.03.02.f"))
clust.in.phy <- subset_samples(clust.in.phy, sample_names(clust.in.phy) != c("myco.06.04.d"))

# all trees with outliers removed
clust.in.bray <- clust.in.phy %>% phyloseq::distance("bray") %>% sqrt
clust.in.nmds <- metaMDS(clust.in.bray, trymax = 200, parallel = 10)
  ## no convergence, run 200 stress 0.207
clust.in.dat <- scores(clust.in.nmds) %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(clust.in.phy) %>% data.frame %>% rownames_to_column("sampID"))
ggplot(clust.in.dat,aes(x=NMDS1,y=NMDS2,fill=Disease)) + #,label = sample_names(clust.in.phy))) +
  geom_point(shape=21,size=3)+
  #geom_text() +
 # stat_ellipse(mapping = aes(color = Disease)) +
  scale_fill_viridis(discrete = T, option = "C") +
  #scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few()

# diseased trees only
clust.dis.bray <- clust.dis.phy %>% phyloseq::distance("bray") %>% sqrt                 
clust.dis.nmds <- metaMDS(clust.dis.bray, trymax = 500, parallel = 10, k=3) 
  ## Converged after 90 runs (with outliers removed), final stress 0.122
clust.dis.dat <-scores(clust.dis.nmds) %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(clust.dis.phy) %>% data.frame %>% rownames_to_column("sampID"))
ggplot(clust.dis.dat,aes(x=NMDS1,y=NMDS2,fill=Disease)) + #, label = sample_names(clust.dis.phy)))+
  geom_point(shape=21,size=3)+
 # stat_ellipse(mapping = aes(color = Disease)) +
  #geom_text() +
  scale_fill_viridis(discrete = T, option = "C") +
  #scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few()
  ## outliers to remove: myco.04.02.e, myco.06.04.d
clust.dis.phy %<>% subset_samples(., sample_names(clust.dis.phy) != c("myco.06.04.d"))
clust.dis.phy %<>% subset_samples(., sample_names(clust.dis.phy) != c("myco.04.02.e"))

health.dis.phy <- subset_samples(clust.prop.phy, Disease==c("Canker", "Healthy"))
health.dis.bray <- health.dis.phy %>% phyloseq::distance("bray") %>% sqrt                 
health.dis.nmds <- metaMDS(health.dis.bray, trymax = 200, parallel = 10) 
  ## no convergence, final stress .2282
health.dis.dat <-scores(health.dis.nmds) %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(health.dis.phy) %>% data.frame %>% rownames_to_column("sampID"))
ggplot(health.dis.dat,aes(x=NMDS1,y=NMDS2,fill=Disease)) + #, label = sample_names(health.dis.phy)))+
  geom_point(shape=21,size=3)+
  stat_ellipse(mapping = aes(color = Disease)) +
  #geom_text() +
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few()
  ## outliers to remove: myco.03.04.a, myco.03.02.f, myco.04.02.e
health.dis.phy %<>% subset_samples(., sample_names(health.dis.phy) != c("myco.03.04.a"))
health.dis.phy %<>% subset_samples(., sample_names(health.dis.phy) != c("myco.03.02.f"))
health.dis.phy %<>% subset_samples(., sample_names(health.dis.phy) != c("myco.04.02.e"))

noncanker.phy <- subset_samples(clustered.phy, Disease==c("Non-canker", "Healthy"))
noncanker.bray <- noncanker.phy %>% phyloseq::distance("bray") %>% sqrt                 
noncanker.nmds <- metaMDS(noncanker.bray, trymax = 200, parallel = 10) 
## no convergence, final stress .198
noncanker.dat <-scores(noncanker.nmds) %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(noncanker.phy) %>% data.frame %>% rownames_to_column("sampID"))
ggplot(noncanker.dat,aes(x=NMDS1,y=NMDS2,fill=Disease, label = sample_names(noncanker.phy)))+ #, label = sample_names(clustered.phy)
  geom_point(shape=21,size=3)+
  geom_text() +
  ggthemes::theme_few()
  ## outliers: myco.07.12.b, myco.08.07.f, myco.05.09.e, myco.06.02.a, 
  ## myco.07.07.h, myco.05.05.d, myco.06.08.a, myco.07.02.d

# trees by source population - subsetting to river sources with 10+ genotypes represented?
pop.bray <- pop.phy %>% phyloseq::distance("bray") %>% sqrt
pop.nmds <- metaMDS(pop.bray, trymax = 200, parallel = 10)
## no convergence, final stress 0.1915
pop.dat <- scores(pop.nmds) %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(pop.phy) %>% data.frame %>% rownames_to_column("sampID"))
ggplot(pop.dat,aes(x=NMDS1,y=NMDS2,fill=Pop))+ #,label = sample_names(pop.phy))) +
  geom_point(shape=21,size=3)+
  #geom_text() +
  #stat_ellipse(mapping = aes(color = Pop)) +
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few()

pop2.5.bray <- pop2.5.phy %>% phyloseq::distance("bray") %>% sqrt
pop2.5.nmds <- metaMDS(pop2.5.bray, trymax = 200, parallel = 10, k=3)
## no convergence, final stress 0.1279
pop2.5.dat <- scores(pop2.5.nmds) %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(pop2.5.phy) %>% data.frame %>% rownames_to_column("sampID"))
ggplot(pop2.5.dat, aes(x=NMDS1, y=NMDS2, shape=Disease, color=Pop)) + #, label = sample_names(pop2.5.phy))) +
  geom_point(size=3)+
  #geom_text() +
  #stat_ellipse(mapping = aes(color = Pop)) +
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few()
  ## outliers: myco.06.04.d, myco.07.05.d, myco.03.02.f
pop2.5.phy <- subset_samples(pop2.5.phy, sample_names(pop2.5.phy) != c("myco.06.04.d")) 
pop2.5.phy <- subset_samples(pop2.5.phy, sample_names(pop2.5.phy) != c("myco.07.05.d")) 
pop2.5.phy <- subset_samples(pop2.5.phy, sample_names(pop2.5.phy) != c("myco.03.02.f")) 

# ordination of core taxa
core.bray <- core.tax.phy %>% phyloseq::distance("bray") %>% sqrt
core.nmds <- metaMDS(core.bray, trymax = 200, parallel = 10, k=3)
  ## did not converge, final stress 0.136
core.nmds.dat <- scores(core.nmds) %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(core.tax.phy) %>% data.frame %>% rownames_to_column("sampID"))
ggplot(core.nmds.dat,aes(x=NMDS1,y=NMDS2,fill=Disease))+ 
  geom_point(shape=21,size=3)+
  stat_ellipse(mapping = aes(color = Disease)) +
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few()

#### PerMANOVA tests ####
# by disease category
clust2.5.dist <- clust2.5.phy %>% phyloseq::distance("bray")  
adonis(clust2.5.dist~sample_data(clust2.5.phy)$Disease)
pair.dis<-pairwise.adonis(clust2.5.dist,factors=sample_data(clust2.5.phy)$Disease)

disease.dist <- clust.dis.phy %>% phyloseq::distance("bray")
adonis(disease.dist~sample_data(clust.dis.phy)$Disease)
# by population
pop2.5.dist <- pop2.5.phy %>% phyloseq::distance("bray")  
adonis(pop2.5.dist~sample_data(pop2.5.phy)$Pop)
# by population, disease and combo
adonis(pop2.5.dist~sample_data(pop2.5.phy)$Disease*sample_data(pop2.5.phy)$Pop)

# core composition by disease category
core.tax.dist <- core.tax.phy %>% phyloseq::distance("bray")
adonis(core.tax.dist~sample_data(core.tax.phy)$Disease)
core.pairs<-pairwise.adonis(core.tax.dist,factors=sample_data(core.tax.phy)$Disease)

#### vector analysis (based on ordinations) ####
  ## vector analysis of all OTUs at 2.5% OTU prevalence, outliers removed
clust.in.otu <- otu_table(clust.in.phy) %>% data.frame()
clust.in.meta <- sample_data(clust.in.phy) %>% data.frame()
clust.in.var <- as.data.frame(scores(clust.in.nmds, display = 'sites'))
clust.in.var <- clust.in.meta %>% cbind(clust.in.var) 
clust.in.var <- clust.in.otu %>% cbind(clust.in.var) 
clust.in.fit <- envfit(clust.in.nmds, clust.in.var, na.rm = TRUE)
#pval.adjust <- p.adjust(clust.fit$pvals, method = "fdr")
#clust.fit$vectors$pvals <- pval.adjust
clust.in.filt <- data.frame(r = clust.in.fit$vectors$r,pvals = clust.in.fit$vectors$pvals) %>%
  rownames_to_column(var = "var") %>%
  filter(r >= 0.2, pvals <= 0.05)

clust.var.scrs <- as.data.frame(scores(clust.in.fit, display = "vectors"))
clust.var.scrs <- cbind(clust.var.scrs, variables = rownames(clust.var.scrs))
clust.var.scrs %<>% filter(variables %in% clust.in.filt$var)

plot <- ggplot(clust.in.dat) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, color = Disease)) +
  coord_fixed() +
  geom_segment(data = clust.var.scrs, size = 2,
               aes(x = 0, xend = NMDS1*.8, y = 0, yend = NMDS2*.8),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  geom_text(data = clust.var.scrs, aes(x = NMDS1, y = NMDS2, label = variables), size = 5, hjust = 1, vjust = 0)

# vector analysis for just OTUs listed in original analysis
  ## OTUs: 1,2,3,8
clust.spp.var <- clust.in.otu %>% dplyr::select("S. musiva" = OTU.1, 
                                             "M. tassiana" = OTU.2, 
                                             "A. pullulans" = OTU.3, 
                                             "C. subcutanea" = OTU.8)
clust.spp.fit <- envfit(clust.in.nmds, clust.spp.var)
clust.spp.filt <- data.frame(r = clust.spp.fit$vectors$r,pvals = clust.spp.fit$vectors$pvals) %>%
  rownames_to_column(var = "var") %>%
  filter(r >= 0.2, pvals <= 0.05)

clust.spp.scrs <- as.data.frame(scores(clust.spp.fit, display = "vectors"))
clust.spp.scrs <- cbind(clust.spp.scrs, variables = rownames(clust.spp.scrs))
clust.spp.scrs %<>% filter(variables %in% clust.spp.filt$var)

spp.var.plot <- ggplot(clust.in.dat) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, color = Disease)) +
  scale_color_viridis(discrete = T, option = "C", "Disease") +
  coord_fixed() +
  geom_segment(data = clust.spp.scrs, size = 2,
               aes(x = 0, xend = NMDS1*.8, y = 0, yend = NMDS2*.8),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  geom_text(data = clust.spp.scrs, aes(x = NMDS1, y = NMDS2, label = variables), 
            size = 3, hjust = 0.5, vjust = 0) +
  ggthemes::theme_few()

#### indicator species analysis ####
disease.cat <- c(rep(1,64), rep(2,135), rep(3,37))
clust.otu <- read.csv("data/clust.otu.csv", as.is = T, row.names = 1)
## note: group 1 is canker, group 2 is healthy, group 3 is non-canker
## note 2: using data.frame of regular OTU, samples ordered by disease category
disease.indval <- multipatt(clust.otu, disease.cat, control = how(nperm = 999))
summary(disease.indval)
