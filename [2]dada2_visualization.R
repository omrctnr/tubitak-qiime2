library(phyloseq)
library(dplyr)
library(tibble)
library(tidyverse)
library(Biostrings)
library(paletteer)
library(microViz)
library(fantaxtic)
library(RColorBrewer) # Special Color Palettes
library(tidyr)
library(ggplot2)
library(ggtext)
library(ggpubr)

# values = getPalette(colourCount)

set.seed(123123)

# load("C:/Users/ASUS/Desktop/prosthesis_microbiota_study/R_Scripts/dada2_visualization_environment.RData")

metadata <- read.table("C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/metadata.txt",
                       header = TRUE, sep = "")

tax_filtered_2 <- tax_filtered[,!names(tax_filtered) %in% c("Sequence")]


metadata <- metadata %>%
  mutate(Condition = replace(Condition, Condition == "AL", "Aseptic Loosening"),
         Condition = replace(Condition, Condition == "PJI", "Prosthetic Joint Infection"),
         dummy_var = seq(1,54),
         Site = replace(Site, Site == "Implant", "Biofilm"))
  

ps <- phyloseq(otu_table(counts_filtered, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(as.matrix(tax_filtered_2)))


# Remove taxa with small mean relative abundance.
# remotes::install_github("vmikk/metagMisc")
library("metagMisc")
ps_filtered <- phyloseq_filter_taxa_rel_abund(ps, frac = 0.001)

sum(ps_filtered@otu_table)/ sum(ps@otu_table) * 100

sum(ps@otu_table) - sum(ps_filtered@otu_table) 

# the rarefaction depth chosen is the minimum sample depth 
# (in this case 604 reads per sample)

ps_rarefied <- rarefy_even_depth(ps_filtered,rngseed=123123, 
                                 sample.size=min(sample_sums(ps_filtered)), replace = FALSE)

## ## ## ## ## ## ##
## VISUALIZATION  ##
## ## ## ## ## ## ##

## general review phyla with stacked bar plot strafied by condition

ps_filt_rel_ab <- transform_sample_counts(ps_filtered, function(x)100* x / sum(x))

general_phylum <- ps_filt_rel_ab %>%
  plot_bar(fill= "Phylum")+
  geom_bar(aes(fill=Phylum), stat="identity", position="stack")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(12,"pt"),
        legend.title = element_text(face = "bold", size = 11),
        strip.text.x = element_text(face = "bold", size = 15),
        title=element_text(face = "bold",size = 13),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(face="bold", colour = "black",size = 12))+
  facet_wrap(~ Condition, scales = "free_x")+
  labs(name="Phylum",x=NULL ,y="   ", title = NULL)+
  guides(x = "none")+
  # scale_fill_manual(name=NULL,values=c(brewer.pal(9,"Paired"),"gray"))
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#33A02C",
                               "#B2DF8A","#E31A1C","#FB9A99",
                               "#FDBF6F","#FF7F00","#CAB2D6"))

general_phylum

ggsave("general_phylum_gg.png",plot = general_phylum, height = 3.5, width = 8.5)


## top 10 family with stacked bar plot strafied by Condition

# select top 10 ASVs in a phyloseq object
top10_family <- top_taxa(ps_filtered, tax_level = "Family", 
                         n_taxa = 10,include_na_taxa = FALSE)

df_top10_family <- top10_family$top_taxa
vec_top10_family <- df_top10_family[["Family"]]

ps_filtered_family <- subset_taxa(ps_filtered, Family %in% vec_top10_family)

ps_filtered_family <- transform_sample_counts(ps_filtered_family, function(x)100* x/sum(x))

top10family<- plot_bar(ps_filtered_family, fill= "Family")+
  geom_bar(aes(), stat="identity", position="stack")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(12,"pt"),
        legend.title =  element_text(size = 11,face = "bold"),
        legend.text = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_markdown(face="bold", colour = "black",size = 15))+
  facet_wrap(~ Condition, scales = "free_x")+
  labs(x=NULL ,y="Relative Abundance <br>(%)", title = NULL)+
  guides(x = "none")+
  scale_fill_manual(name="Family",values=c(brewer.pal(10,"Paired"),"gray"))

top10family

ggsave("top10family_gg.png",plot = top10family, height = 3.5, width = 8.5)

###

## top 10 genus with stacked bar plot strafied by Condition

# select top 10 ASVs in a phyloseq object

top15_Genus <- top_taxa(ps_filtered, tax_level = "Genus", 
                        n_taxa = 15,include_na_taxa = FALSE)

df_top15_Genus <- top15_Genus$top_taxa
vec_top15_Genus <- df_top15_Genus[["Genus"]]

ps_filtered_Genus <- subset_taxa(ps_filtered, Genus %in% vec_top15_Genus)

ps_filtered_Genus <- transform_sample_counts(ps_filtered_Genus, function(x)100* x/sum(x))

top15Genus<- plot_bar(ps_filtered_Genus, fill= "Genus")+
  geom_bar(aes(), stat="identity", position="stack")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(12,"pt"),
        legend.title =  element_text(size = 11,face = "bold"),
        legend.text = element_text(face = "italic"),
        title=element_text(face = "bold",size = 13),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text(face="bold", colour = "black",size = 12))+
  facet_wrap(~ Condition, scales = "free_x")+
  labs(x=NULL ,y="   ", title = NULL)+
  guides(x = "none")+
  scale_fill_manual(name="Genus",values = rev(distinct_palette(n = 16, pal = "kelly",add = NA)))

top15Genus

ggsave("top15Genus_gg.png",plot = top15Genus, height = 3.5, width = 8.5)

##### merge plots
library(cowplot) # merging function

general_phylum
top10family
top15Genus

d1<- plot_grid(general_phylum,top10family,top15Genus, nrow = 3,ncol = 1,labels = c("A","B","C"),
               label_size = 15)


ggsave("merged_test.tiff", plot = d1, height=8.5, width=8.5)

####### General review box plots

# phylum

ps_filtered_phylum <- tax_glom(ps_filt_rel_ab, "Phylum")
taxa_names(ps_filtered_phylum) <- tax_table(ps_filtered_phylum)[, "Phylum"]

phylum_review <- psmelt(ps_filtered_phylum)%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  theme_bw()+
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, scales = "free")+
  theme(legend.position = "none")

## testing the relative abundances and adding the significance 

ps_phylum_rev_melt_df <- psmelt(ps_filtered_phylum)%>%
  mutate(Condition = replace(Condition, Condition == "Aseptic Loosening", "AL"),
         Condition = replace(Condition, Condition == "Prosthetic Joint Infection","PJI"))

# testing rel abundance 
phylum.rel.stat <- ps_phylum_rev_melt_df  %>% 
  group_by(OTU) %>% 
  drop_na()%>%
  wilcox_test(data = ., Abundance ~ Condition) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")%>%
  add_xy_position()%>%
  mutate(y.position= y.position - 13)

phylum_review <- phylum_review + stat_pvalue_manual(phylum.rel.stat, label = "p.adj.signif")

ggsave("phylum_review_gg.png",plot = phylum_review, height = 8, width = 10)

# family

top15_family <- top_taxa(ps_filtered, tax_level = "Family", n_taxa = 10,include_na_taxa = FALSE)
df_top15_family <- top15_family$top_taxa
vec_top15_family <- df_top15_family[["Family"]]

ps_filt_rel_ab <- transform_sample_counts(ps_filtered, function(x)100* x/sum(x))

ps_filtered_family <- tax_glom(ps_filt_rel_ab, "Family")
taxa_names(ps_filtered_family) <- tax_table(ps_filtered_family)[, "Family"]

ps_filtered_family <- subset_taxa(ps_filtered_family, Family %in% vec_top15_family)

family_review <- psmelt(ps_filtered_family)%>%
  mutate(Condition = replace(Condition, Condition == "Aseptic Loosening", "AL"),
         Condition = replace(Condition, Condition == "Prosthetic Joint Infection","PJI"))%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  theme_bw()+
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, scales = "free", nrow = 3,ncol = 5)+
  theme(legend.position = "none")

## testing the relative abundances and adding the significance 

ps_family_rev_melt_df <- psmelt(ps_filtered_family)%>%
  mutate(Condition = replace(Condition, Condition == "Aseptic Loosening", "AL"),
         Condition = replace(Condition, Condition == "Prosthetic Joint Infection","PJI"))

# testing rel abundance 
family.rel.stat <- ps_family_rev_melt_df  %>% 
  group_by(OTU) %>% 
  drop_na()%>%
  wilcox_test(data = ., Abundance ~ Condition) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")%>%
  add_xy_position()%>%
  mutate(y.position= y.position - 13)

family_review <- family_review + stat_pvalue_manual(family.rel.stat, label = "p.adj.signif")

ggsave("family_review_gg.png",plot = family_review, height = 8, width = 10)


# genus

top15_Genus <- top_taxa(ps_filtered, tax_level = "Genus", n_taxa = 15,include_na_taxa = FALSE)
df_top15_Genus <- top15_Genus$top_taxa
vec_top15_Genus <- df_top15_Genus[["Genus"]]

ps_filt_rel_ab <- transform_sample_counts(ps_filtered, function(x)100* x/sum(x))

ps_filtered_Genus <- tax_glom(ps_filt_rel_ab, "Genus")
taxa_names(ps_filtered_Genus) <- tax_table(ps_filtered_Genus)[, "Genus"]

ps_filtered_Genus <- subset_taxa(ps_filtered_Genus, Genus %in% vec_top15_Genus)

Genus_review <- psmelt(ps_filtered_Genus)%>%
  mutate(Condition = replace(Condition, Condition == "Aseptic Loosening", "AL"),
         Condition = replace(Condition, Condition == "Prosthetic Joint Infection","PJI"))%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  theme_bw()+
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, scales = "free", nrow = 3,ncol = 5)+
  theme(legend.position = "none")

## testing the relative abundances and adding the significance 

ps_gen_rev_melt_df <- psmelt(ps_filtered_Genus)%>%
  mutate(Condition = replace(Condition, Condition == "Aseptic Loosening", "AL"),
         Condition = replace(Condition, Condition == "Prosthetic Joint Infection","PJI"))

# testing rel abundance 
gen.rel.stat <- ps_gen_rev_melt_df  %>% 
                  group_by(OTU) %>%
                  drop_na()%>%
                  wilcox_test(data = ., Abundance ~ Condition) %>% 
                  adjust_pvalue(method = "fdr") %>%
                  add_significance("p.adj")%>%
                  add_xy_position()%>%
                  mutate(y.position= y.position - 13)

Genus_review <- Genus_review + stat_pvalue_manual(gen.rel.stat, label = "p.adj.signif")

ggsave("Genus_review_gg.png",plot = Genus_review, height = 8, width = 10)
