library(tibble)
library(Biostrings)
library(microbiomeMarker) # differential abundance analysis
library(tidyr)
library(stringr)
library(ggplot2)
library(phyloseq)
library(dplyr) # data manipulation
library(tidyverse) # data manipulation
library(rstatix) # statistical test
library(ggpubr) # color palette and theme settings
library(cowplot) # merges ggplot objects
library(vegan) # color palette
library(ape)    # pcoa() function
library(ggtext) # provides the writting options based on html | italic | element_markdown()
library(ggsci) # color palette and theme settings
library(pairwiseAdonis) # statistical test
library(RColorBrewer) # color palette
library(relayer) # extra scale in ggplot2 >> rename_geom_aes() function
library(ggnewscale) # extra scale in ggplot2
library(Polychrome) # color palette
library(ggthemes) # color palette
library(paletteer) # color palette
library(microViz)

set.seed(16022000)

ps<-readRDS("PJI_forward_ps_object2.rds")

# Reorder groups and samples
ps@sam_data$Condition <- factor(ps@sam_data$Condition,levels = (c("AL","PJI")))
ps@sam_data$Site[ps@sam_data$Site =="Implant"] <- "Biofilm"
ps@sam_data$Site[ps@sam_data$Site =="Tissue"] <- "Deep Tissue"
ps@sam_data$Site[ps@sam_data$Site =="Synovial_Fluid"] <- "Synovial Fluid"
ps@sam_data$Site <- factor(ps@sam_data$Site,
                                           levels = (c("Biofilm","Deep Tissue","Synovial Fluid")))

#check sample read sums
as.data.frame(sort(sample_sums(ps),decreasing = TRUE))

# filter out taxa with mean relative abundance less than 0.0001 and
ps_relab <- transform_sample_counts(ps, function(x)100* x / sum(x))
ps_rm_taxa<-subset_taxa(ps_relab, !(Order %in% "Chloroplast" | Family %in% "Mitochondria" | Phylum %in% c(NA)| Kingdom %in% c("Eukaryota","Archaea",NA)))
ps_rm_taxa<-filter_taxa(ps_rm_taxa, function(x) mean(x) > 0.001*100, TRUE) 

# create new phyloseq object according to the removed taxa
myTaxa<-rownames(ps_rm_taxa@tax_table)
ps_filtered<-prune_taxa(myTaxa, ps)
ps_filtered

# check samples and see how many reads counts were filtered out
as.data.frame(sort(sample_sums(ps_filtered),decreasing = TRUE))
sum(ps_filtered@otu_table) / sum(ps@otu_table)

# rarefy ps_filtered object in order to calculate alpha and beta diversity
ps_rarefied<-rarefy_even_depth(ps_filtered, rngseed=16022000, 
                               sample.size=min(sample_sums(ps_filtered)))

# remove redundant variables
rm(ps_relab, myTaxa)
###### Tax bar plots ####
## top 15 genus with stacked bar plot strafied by Condition
new_tax_table<-as.data.frame(tax_table(ps_filtered))
new_tax_table$Genus[new_tax_table$Genus =="Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia-<br>Caballeronia-Paraburkholderia"

tax_table(ps_filtered) <- as.matrix(new_tax_table)

ps_genus_glom<-tax_glom(ps_filtered, "Genus")

top10_genus_asv<-names(sort(taxa_sums(ps_genus_glom)[1:15],decreasing = TRUE))

top10_genus_taxa<-as.data.frame(ps_genus_glom@tax_table[top10_genus_asv,])[["Genus"]]

ps_filtered_genus <- subset_taxa(ps_filtered, Genus %in% top10_genus_taxa)

ps_filtered_genus <- transform_sample_counts(ps_filtered_genus, function(x)100* x/sum(x))

top15Genus_color <- set_names(top10_genus_taxa,)

top15Genus<- plot_bar(ps_filtered_genus, fill= "Genus")+
  geom_bar(aes(), stat="identity", position="stack",colour="White")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(10,"pt"),
        legend.title = element_text(face = "bold", size = 10),
        strip.text.x = element_text(face = "bold", size = 11),
        title=element_text(face = "bold",size = 11),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_text(face="bold", colour = "black",size = 7),
        axis.ticks.x= ,
        axis.title.y = element_text(face="bold", colour = "black",size = 10),
        legend.text = element_markdown(face = "italic"),
        legend.position = "none")+
  facet_wrap(~ Condition, scales = "free_x")+
  labs(x="Sample" ,y="Relative Abundance (%)", title = NULL)+
  guides(x = "none")+
  scale_fill_manual(name="Genus",values = rev(distinct_palette(n = 16, pal = "kelly",add = NA)))

top15Genus

# ggsave("top15Genus_gg.png",plot = top15Genus, height = 3.5, width = 8.5)

## merged genus

ps_filtered_merged<-merge_samples(ps_filtered,group = "Condition")
ps_filtered_merged_genus <- subset_taxa(ps_filtered_merged, Genus %in% top10_genus_taxa)
ps_filtered_merged_genus_relab <- transform_sample_counts(ps_filtered_merged_genus, function(x)100* x / sum(x))


top15Genus_merged<- plot_bar(ps_filtered_merged_genus_relab, fill= "Genus")+
  geom_bar(aes(), stat="identity", position="stack",colour="White")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(10,"pt"),
        legend.title = element_text(face = "bold", size = 10),
        strip.text.x = element_text(face = "bold", size = 11),
        title=element_text(face = "bold",size = 11),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_text(face="bold", colour = "black",size = 10),
        axis.ticks.x= ,
        axis.title.y = element_text(face="bold", colour = "black",size = 10),
        legend.text = element_markdown(face = "italic"),
        legend.position = "none")+
  # facet_wrap(~ Condition, scales = "free_x")+
  labs(x=NULL ,y="Relative Abundance (%)", title = NULL)+
  # guides(x = "none")+
  scale_fill_manual(name="Genus",values = rev(distinct_palette(n = 16, pal = "kelly",add = NA)))


top15Genus_merged
# ggsave("top15Genus_merged.png",plot = top15Genus_merged, height = 3.5, width = 5.5)

### merged genera v2
ps_filtered_merged<-merge_samples(ps_filtered,group = "Condition")
ps_filtered_merged_genus <- subset_taxa(ps_filtered_merged, Genus %in% top10_genus_taxa)
ps_filtered_merged_genus_relab <- transform_sample_counts(ps_filtered_merged_genus, function(x)100* x / sum(x))


top15Genus_mergedv2<- ps_filtered%>%
  subset_taxa(Genus %in% top10_genus_taxa)%>%
  psmelt()%>%
  as.data.frame()%>%
  group_by(Site,Condition)%>%
  mutate(Abundance = Abundance / sum(Abundance)*100)%>%
  # mutate(Genus= factor(Genus,levels = top10_genus_taxa))%>%
  
  ggplot(aes(fill=Genus,x=Site,y=Abundance))+
  geom_bar(aes(), stat="identity", position="stack")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(10,"pt"),
        legend.title = element_text(face = "bold", size = 10),
        strip.text.x = element_text(face = "bold", size = 11),
        title=element_text(face = "bold",size = 11),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_text(face="bold", colour = "black",size = 7),
        axis.ticks.x= ,
        axis.title.y = element_text(face="bold", colour = "black",size = 10),
        legend.text = element_markdown(face = "italic"),
        legend.position = "none")+
  facet_wrap(~ Condition, scales = "free_x")+
  labs(x=NULL ,y="Relative Abundance (%)", title = NULL)+
  # guides(x = "none")+
  scale_fill_manual(name="Genus",values = rev(distinct_palette(n = 16, pal = "kelly",add = NA)))


top15Genus_mergedv2
# ggsave("top15Genus_mergedv2.png",plot = top15Genus_mergedv2, height = 3.5, width = 7.5)

### merge plots###

top15Genus_legend<- plot_bar(ps_filtered_genus, fill= "Genus")+
  geom_bar(aes(), stat="identity", position="stack",colour="White")+geom_col(color=NA)+
  theme(
    legend.text = element_markdown(size = 10,face = "italic"),
    legend.key.size = unit(12,"pt"),
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.position = "top",
    panel.border = element_rect(colour = "black"))+
  scale_fill_manual(name="Genus",values = rev(distinct_palette(n = 16, pal = "kelly",add = NA)))+
  guides(fill = guide_legend(title.position = "top"))

as_ggplot(get_legend(top15Genus_legend))

merged_genera <- plot_grid(get_legend(top15Genus_legend),
                           top15Genus,labels = c("","A"),
                           plot_grid(top15Genus_merged,top15Genus_mergedv2,labels = c("B","C"),rel_widths = c(1,1.5)),
                           nrow = 3,ncol = 1,
                           rel_heights = c(1,2.5,2.5))+
  theme(plot.background = element_rect(fill = "White"))

ggsave(plot = merged_genera, "merged_genera.tiff",height = 8 ,width =8)


##### beta diversity #####
set.seed(16022000)

counts_df <- as.data.frame(ps_rarefied@otu_table)
metadata <- data.frame(ps_rarefied@sam_data)

# create Bray-Curtis Distance
bc_distance <- vegdist(counts_df, method = "bray")
bc_pcoa <- pcoa(bc_distance)

# get principal coordinate values
metadata$pcoa1 <- bc_pcoa$vectors[,1]
metadata$pcoa2 <- bc_pcoa$vectors[,2]
metadata$pcoa3 <- bc_pcoa$vectors[,3]

# get variances of principal coordinate values 
bc_pcoa$values$Relative_eig[[1]] * 100
bc_pcoa$values$Relative_eig[[2]] * 100

# PERMANOVA Test according to the day and was significant
permanova_cond <- adonis2(bc_distance ~ Condition, data=metadata,permutations = 9999)
permanova_cond

permanova_site <- adonis2(bc_distance ~ Site, data=metadata,permutations = 9999)
permanova_site

condition_pcoa_plot <- metadata %>% 
  ggplot(aes(x=pcoa1, y=pcoa2, color= Condition), show.legend=FALSE)+
  geom_point(size = 2) + 
  stat_ellipse(linetype=1,lwd=0.75) +
  theme_cowplot()+
  theme(axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        panel.background= element_rect(colour = "black"))+
  labs(x = paste("PC1",paste0("[",round(bc_pcoa$values$Relative_eig[[1]] * 100,digits = 2),"]")),
       y = paste("PC2",paste0("[",round(bc_pcoa$values$Relative_eig[[2]] * 100,digits = 2),"]")))+
  scale_color_lancet(name="Diagnosis")
condition_pcoa_plot

# ggsave("condition_pcoa_plot.png",plot=condition_pcoa_plot, height=4, width=6)  


site_pcoa_plot <- metadata %>% 
  ggplot(aes(x=pcoa1, y=pcoa2, color= Site), show.legend=FALSE)+
  geom_point(size = 2) + 
  stat_ellipse(linetype=1,lwd=0.75) +
  theme_cowplot()+
  theme(axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        panel.background= element_rect(colour = "black"))+
  labs(x = paste("PC1",paste0("[",round(bc_pcoa$values$Relative_eig[[1]] * 100,digits = 2),"]")),
       y = paste("PC2",paste0("[",round(bc_pcoa$values$Relative_eig[[2]] * 100,digits = 2),"]")))+
  scale_color_lancet(name="Condition")
site_pcoa_plot
# ggsave("site_pcoa_plot.png",plot=site_pcoa_plot, height=4, width=6) 


combined_pcoa_plot <- metadata %>% 
  ggplot(aes(x=pcoa1, y=pcoa2, color= Condition, shape=Site), show.legend=FALSE)+
  geom_point(aes(),size = 2) +
  theme_cowplot()+
  theme(axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        panel.background= element_rect(colour = "black"))+
  labs(x = paste("PC1",paste0("[",round(bc_pcoa$values$Relative_eig[[1]] * 100,digits = 2),"]")),
       y = paste("PC2",paste0("[",round(bc_pcoa$values$Relative_eig[[2]] * 100,digits = 2),"]")))+
  scale_shape_manual(values=c(15,17,3),name="Site")+
  scale_color_lancet(name="Diagnosis")
combined_pcoa_plot

# ggsave("Combined_PCoA_plot.tiff",plot=combined_pcoa_plot, height=4, width=6)

###### Differential Abundance #####
# https://ycl6.github.io/16S-Demo/index.html
# https://www.sciencedirect.com/science/article/pii/S2352304217300351
set.seed(16022000)
lefse_diff_abd <-run_lefse(ps_rarefied,
                           wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, kw_cutoff = 0.05,lda_cutoff = 1,
                           group = 'Condition')

my_Df_condition<-as.data.frame(lefse_diff_abd@marker_table)
my_Df_condition<-my_Df_condition%>%
  data.frame()%>%
  separate(feature, sep = "\\|", remove = FALSE, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"))%>%
  mutate(feature_main = case_when(str_detect(feature, "s__") ~ str_replace_all(str_extract(feature, "s__.*"), "__", "_"),
                                  !str_detect(feature, "s__") & str_detect(feature, "g__")~ str_replace_all(str_extract(feature, "g__.*"), "__", "_"),
                                  !str_detect(feature, "g__") & str_detect(feature, "f__")~ str_replace_all(str_extract(feature, "f__.*"), "__", "_"),
                                  !str_detect(feature, "f__") & str_detect(feature, "o__")~ str_replace_all(str_extract(feature, "o__.*"), "__", "_"),
                                  !str_detect(feature, "o__") & str_detect(feature, "c__")~ str_replace_all(str_extract(feature, "c__.*"), "__", "_"),
                                  !str_detect(feature, "c__") & str_detect(feature, "p__")~ str_replace_all(str_extract(feature, "p__.*"), "__", "_"),
                                  !str_detect(feature, "p__") & str_detect(feature, "k__")~ str_replace_all(str_extract(feature, "k__.*"), "__", "_"),
                                  TRUE ~ NA_character_))

my_Df_condition$feature_main <- factor(my_Df_condition$feature_main, levels = my_Df_condition$feature_main)

remove_tx <- c("s_Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium_s_",
               "g_Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
               "o_Veillonellales-Selenomonadales",
               "o_Lactobacillales",
               "o_Pasteurellales")

my_Df_condition <-filter(my_Df_condition, !(feature_main %in% remove_tx))


lefse_plot<-my_Df_condition%>%
  mutate(feature_main = sub("(.{2})(.*)", "\\1*\\2", feature_main),
         feature_main = paste(feature_main,"*", sep = ""))%>%
  filter(ef_lda>4)%>%
  mutate(ef_lda = if_else(enrich_group=="AL", -1 * ef_lda, ef_lda),
         feature_main = fct_reorder(feature_main, ef_lda))%>%
  ggplot(aes(x=feature_main,y=ef_lda))+
  geom_bar(aes(fill=enrich_group),stat = 'identity', col="black",width = 0.6)+
  scale_fill_manual(values = rev(ggsci::pal_nejm(palette = "default")(2)[1:2]) ,name="Diagnosis",)+
  labs(y="LDA Score (log 10)", x=NULL, title = NULL)+
  theme_minimal() +
  coord_flip()+
  scale_y_continuous(expression(log[10](italic("LDA Score"))),
                     breaks = seq(-5,5, by=2), limits = c(-5, 5))+
  geom_richtext(aes(y = 0,label = feature_main, hjust = ifelse(ef_lda < 0, -0.03, 1.03)),
                size=4,fill = NA, label.colour = NA)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face="bold"),
        legend.key.size = unit(15,"pt"),
        legend.position = "top",
        legend.justification = 0.05,
        legend.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank(),
        plot.background = element_rect(),
        plot.caption = element_text(face = "italic", size = 12),
  )

lefse_plot

# ggsave(plot = lefse_plot,filename = "lefse_condition.jpg",height = 5,width = 6)

#lefse site

ps_trans<-subset_samples(ps_rarefied,Site %in% c("Biofilm","Deep Tissue"))
ps_rareifed_pjı_bio_deep<-subset_samples(ps_trans,Condition %in% c("PJI"))

set.seed(16022000)
lefse_diff_abd_site <-run_lefse(ps_rareifed_pjı_bio_deep,
                                wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, kw_cutoff = 0.05,lda_cutoff = 1,
                                group = 'Site')

my_Df_Site<-as.data.frame(lefse_diff_abd_site@marker_table)
my_Df_Site<-my_Df_Site%>%
  data.frame()%>%
  separate(feature, sep = "\\|", remove = FALSE, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"))%>%
  mutate(feature_main = case_when(str_detect(feature, "s__") ~ str_replace_all(str_extract(feature, "s__.*"), "__", "_"),
                                  !str_detect(feature, "s__") & str_detect(feature, "g__")~ str_replace_all(str_extract(feature, "g__.*"), "__", "_"),
                                  !str_detect(feature, "g__") & str_detect(feature, "f__")~ str_replace_all(str_extract(feature, "f__.*"), "__", "_"),
                                  !str_detect(feature, "f__") & str_detect(feature, "o__")~ str_replace_all(str_extract(feature, "o__.*"), "__", "_"),
                                  !str_detect(feature, "o__") & str_detect(feature, "c__")~ str_replace_all(str_extract(feature, "c__.*"), "__", "_"),
                                  !str_detect(feature, "c__") & str_detect(feature, "p__")~ str_replace_all(str_extract(feature, "p__.*"), "__", "_"),
                                  !str_detect(feature, "p__") & str_detect(feature, "k__")~ str_replace_all(str_extract(feature, "k__.*"), "__", "_"),
                                  TRUE ~ NA_character_))

my_Df_Site$feature_main <- factor(my_Df_Site$feature_main, levels = my_Df_Site$feature_main)

remove_tx2 <- c("o_Corynebacteriales")

my_Df_Site <-filter(my_Df_Site, !(feature_main %in% remove_tx2))

lefse_plot_site<-my_Df_Site%>%
  mutate(feature_main = sub("(.{2})(.*)", "\\1*\\2", feature_main),
         feature_main = paste(feature_main,"*", sep = ""))%>%
  filter(ef_lda>4)%>%
  mutate(ef_lda = if_else(enrich_group=="Deep Tissue", -1 * ef_lda, ef_lda),
         feature_main = fct_reorder(feature_main, ef_lda))%>%
  ggplot(aes(x=feature_main,y=ef_lda))+
  geom_bar(aes(fill=enrich_group),stat = 'identity', col="black",width = 0.6)+
  scale_fill_manual(values = ggsci::pal_nejm(palette = "default")(4)[3:4] ,name="Site*",)+
  labs(y="LDA Score (log 10)", x=NULL, title = NULL,caption = "*Only PJI samples were involved in the analysis,<br>excluding synovial fluid specimens")+
  theme_minimal() +
  coord_flip()+
  scale_y_continuous(expression(log[10](italic("LDA Score"))),
                     breaks = seq(-5,5, by=2), limits = c(-5, 5))+
  geom_richtext(aes(y = 0,label = feature_main, hjust = ifelse(ef_lda < 0, -0.03, 1.03)),
                size=4,fill = NA, label.colour = NA)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face="bold"),
        legend.key.size = unit(15,"pt"),
        legend.position = "top",
        legend.justification = 0.05,
        legend.title = element_text(face="bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank(),
        plot.background = element_rect(),
        plot.caption = element_markdown(face = "italic", size = 10)
  )

lefse_plot_site

# ggsave(plot = lefse_plot_site,filename = "lefse_site.jpg",height = 5,width = 6)

# merge plots
d1<- plot_grid(lefse_plot,lefse_plot_site, nrow = 1,ncol = 2,labels = c("A","B"),
               label_size = 13, rel_widths = c(1.2,1))

ggsave("lefse_plots_merged.tiff",plot = d1,height = 5, width = 10)

##### signifcant rel ab box plots ####
species_rel_test<- ps_rm_taxa  %>% 
  psmelt() %>%
  group_by(Species) %>%
  wilcox_test(data = ., Abundance ~ Condition) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")%>%
  add_xy_position(scales = "free")%>%
  filter(p.adj < 0.05)

specie_name<-c('gordonii' = "Streptococcus gordonii")

species_boxplot<-psmelt(ps_rm_taxa)%>%
  filter(Species == "gordonii")%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot() +
  # geom_jitter(aes(color = Species), height = 0, width = .2) +
  theme_bw()+
  labs(x = "", y = "") +
  facet_wrap(~ Species, scales = "free", nrow = 3,ncol = 5,labeller = as_labeller(specie_name))+
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"))+
  stat_pvalue_manual(species_rel_test)

species_boxplot

#2
ps_genus_glom<-tax_glom(ps_filtered, "Genus")

top10_genus_asv<-names(sort(taxa_sums(ps_genus_glom)[1:15],decreasing = TRUE))

top10_genus_taxa<-as.data.frame(ps_genus_glom@tax_table[top10_genus_asv,])[["Genus"]]

genus_rel_test <- ps_rm_taxa  %>% 
  psmelt() %>%
  group_by(Genus) %>%
  wilcox_test(data = ., Abundance ~ Condition) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")%>%
  add_xy_position(scales = "free")%>%
  filter(Genus %in% top10_genus_taxa)

genus_specific1<-psmelt(ps_rm_taxa)%>%
  filter(Genus == filter(genus_rel_test,p.adj<0.05)$Genus[1] )%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot() +
  theme_bw()+
  labs(x = "", y = "") +
  facet_wrap(~ Genus, scales = "free", nrow = 3,ncol = 5)+
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"))+
  stat_pvalue_manual(filter(genus_rel_test,p.adj<0.05)[1,])
#3
genus_specific2<-psmelt(ps_rm_taxa)%>%
  filter(Genus == filter(genus_rel_test,p.adj<0.05)$Genus[2] )%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot() +
  theme_bw()+
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ Genus, scales = "free", nrow = 3,ncol = 5)+
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"))+
  stat_pvalue_manual(filter(genus_rel_test,p.adj<0.05)[2,])

#4

ps_family_glom<-tax_glom(ps_filtered, "Family")

top10_family_asv<-names(sort(taxa_sums(ps_family_glom)[1:10],decreasing = TRUE))

top10_family_taxa<-as.data.frame(ps_family_glom@tax_table[top10_family_asv,])[["Family"]]

family_rel_test <- ps_rm_taxa  %>% 
  psmelt()%>%
  group_by(Family) %>% 
  wilcox_test(data = ., Abundance ~ Condition) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")%>%
  add_xy_position(scales = "free")%>%
  filter(Family %in% top10_family_taxa)

family_specific<-psmelt(ps_rm_taxa)%>%
  filter(Family == filter(family_rel_test,p.adj<0.05)$Family )%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot() +
  theme_bw()+
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ Family, scales = "free", nrow = 3,ncol = 5)+
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"))+
  stat_pvalue_manual(filter(family_rel_test,p.adj<0.05))

merged_signif_taxa_plot<-plot_grid(family_specific,genus_specific1, genus_specific2,species_boxplot)

# ggsave(plot = merged_signif_taxa_plot, "merged_signif_taxa_plot.tiff",height = 5,width = 6)

# combined figure 2-3-4
combined_pcoa_plot
merged_signif_taxa_plot

d1<- plot_grid(lefse_plot,lefse_plot_site, nrow = 1,ncol = 2,labels = c("C1","C2"),
               label_size = 12, rel_widths = c(1.2,1),label_x = 0.01)



combined_f2_f3_f4<-plot_grid(plot_grid(merged_signif_taxa_plot,combined_pcoa_plot,rel_widths = c(1,1.25),
                                       labels = c("A","B"),label_size = 15,label_x = -0.01),d1,ncol = 1,nrow = 2,scale = 0.95,
                             rel_heights = c(1.1,1),labels = c(" ","C"),label_size = 15,label_x = -0.004)+
  theme(plot.background = element_rect(fill = "White"))

ggsave(plot = combined_f2_f3_f4, "combined_f2_f3_f4.tiff",height = 10,width = 10)