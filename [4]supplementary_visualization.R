### stacked bar plots #######
ps_filt_relab <- transform_sample_counts(ps_filtered, function(x)100* x / sum(x))

general_phylum <- ps_filt_relab %>%
  plot_bar(fill= "Phylum")+
  geom_bar(aes(fill=Phylum), stat="identity", position="stack",colour="White")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(10,"pt"),
        legend.title = element_text(face = "bold", size = 11),
        strip.text.x = element_text(face = "bold", size = 11),
        title=element_text(face = "bold",size = 11),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(face="bold", colour = "black",size = 10),
        legend.text = element_text(face = "italic"))+
  facet_wrap(~ Condition, scales = "free_x")+
  labs(name="Phylum",x="" ,y="Relative Abundance (%)", title = NULL)+
  guides(x = "none")+
  scale_fill_brewer(palette = "Paired")
# # scale_fill_manual(name=NULL,values=c(brewer.pal(9,"Paired"),"gray"))
# scale_fill_manual(values = c("#A6CEE3","#1F78B4","#33A02C",
#                              "#B2DF8A","#E31A1C","#FB9A99",
#                              "#FDBF6F","#FF7F00","#CAB2D6"))
general_phylum

# ggsave("general_phylum_gg.tiff",plot = general_phylum, height = 3.5, width = 8.5)

## merged phyla v1
ps_filtered_merged<-merge_samples(ps_filtered,group = "Condition")
ps_filtered_merged_relab <- transform_sample_counts(ps_filtered_merged, function(x)100* x / sum(x))

general_phylum_merged<-ps_filtered_merged_relab %>%
  plot_bar(fill= "Phylum")+
  geom_bar(aes(fill=Phylum), stat="identity", position="stack",colour="White")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(10,"pt"),
        legend.title = element_text(face = "bold", size = 10),
        strip.text.x = element_text(face = "bold", size = 11),
        title=element_text(face = "bold",size = 11),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_text(face="bold", colour = "black",size = 10),
        axis.ticks.x= ,
        axis.title.y = element_text(face="bold", colour = "black",size = 10),
        legend.text = element_text(face = "italic"),
        legend.position = "none",
        plot.title = element_text(size = 12,face = "bold",hjust = 0.5))+
  labs(name="Phylum",x=NULL ,y="Relative Abundance (%)", title = "Phylum")+
  # guides(x = "none")+
  scale_fill_brewer(palette = "Paired")

general_phylum_merged

# ggsave("general_phylum_merged.png",plot = general_phylum_merged, height = 3.5, width = 4.5)

## merged phyla v2
ps_phylum_glom<-tax_glom(ps_filtered, "Phylum")

top10_phylum_asv<-names(sort(taxa_sums(ps_phylum_glom)[1:11],decreasing = TRUE))
top10_phylum_taxa<-as.data.frame(ps_phylum_glom@tax_table[top10_phylum_asv,])[["Phylum"]]

general_phylum_mergedv2<-ps_filtered %>%
  psmelt()%>%
  as.data.frame()%>%
  group_by(Site,Condition)%>%
  mutate(Abundance = Abundance / sum(Abundance)*100)%>%
  # mutate(Phylum= factor(Phylum,levels = top10_phylum_taxa ))%>%
  ggplot(aes(fill=Phylum,x=Site,y=Abundance))+
  geom_bar(aes(), stat="identity", position="stack")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(12,"pt"),
        legend.title = element_text(face = "bold", size = 11),
        strip.text.x = element_text(face = "bold", size = 13),
        title=element_text(face = "bold",size = 13),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_text(face="bold", colour = "black",size = 8),
        # axis.ticks.x=element_blank(),
        axis.title.y = element_text(face="bold", colour = "black",size = 12),
        legend.text = element_text(face = "italic"),
        legend.position = "none",
        plot.title = element_text(size = 12,face = "bold",hjust = 0.5))+
  facet_wrap(~Condition, scales = "free_x",)+
  labs(name="Phylum",x=NULL ,y="Relative Abundance (%)", title = "Phylum")+
  # guides(x = "none")+
  scale_fill_brewer(palette = "Paired")

general_phylum_mergedv2

# ggsave("general_phylum_mergedv2.png",plot = general_phylum_mergedv2, height = 3.5, width = 7.5)

## top 10 family with stacked bar plot strafied by Condition

# select top 10 ASVs in a phyloseq object
ps_family_glom<-tax_glom(ps_filtered, "Family")

top10_family_asv<-names(sort(taxa_sums(ps_family_glom)[1:10],decreasing = TRUE))

top10_family_taxa<-as.data.frame(ps_family_glom@tax_table[top10_family_asv,])[["Family"]]


ps_filtered_family <- subset_taxa(ps_filtered, Family %in% top10_family_taxa)

ps_filtered_family <- transform_sample_counts(ps_filtered_family, function(x)100* x/sum(x))

top10family<- plot_bar(ps_filtered_family, fill= "Family")+
  geom_bar(aes(), stat="identity", position="stack",colour="White")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(10,"pt"),
        legend.title = element_text(face = "bold", size = 11),
        strip.text.x = element_blank(),
        title=element_text(face = "bold",size = 11),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(face="bold", colour = "black",size = 10),
        legend.text = element_text(face = "italic"),
        strip.background = element_blank())+
  facet_wrap(~ Condition, scales = "free_x")+
  labs(x=NULL ,y="Relative Abundance (%)", title = NULL)+
  guides(x = "none")+
  scale_fill_manual(name="Family",values=ggsci::pal_d3(palette = "category10")(10))

top10family

# ggsave("top10family_gg.png",plot = top10family, height = 3.5, width = 8.5)



### merged family

ps_filtered_merged<-merge_samples(ps_filtered,group = "Condition")
ps_filtered_merged_family <- subset_taxa(ps_filtered_merged, Family %in% top10_family_taxa)
ps_filtered_merged_family_relab <- transform_sample_counts(ps_filtered_merged_family, function(x)100* x / sum(x))

top10family_merged<-plot_bar(ps_filtered_merged_family_relab, fill= "Family")+
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
        legend.text = element_text(face = "italic"),
        legend.position = "none",
        plot.title = element_text(size = 12,face = "bold",hjust = 0.5))+
  # facet_wrap(~ Condition, scales = "free_x")+
  labs(x=NULL ,y="Relative Abundance (%)", title = "Family")+
  # guides(x = "none")+
  scale_fill_manual(name="Family",values=ggsci::pal_d3(palette = "category10")(10))

top10family_merged

# ggsave("top10family_merged.png",plot = top10family_merged, height = 3.5, width = 4.5)

## merged family v2

# physeqPhylum = tax_glom(ps_filtered, "Family")
# physeqPhylumRA = transform_sample_counts(physeqPhylum, function(x) x/sum(x))
# physeqPhylumRAF = filter_taxa(physeqPhylumRA, function(x) mean(x) > 0.01, TRUE)
# # Define the vector of phyla names that are still present after your filter
# keepPhyla = get_taxa_unique(physeqPhylumRAF, "Family")
# # Like in the first question, subset to just the phyla that you want, using the original phyloseq object
# physeqF = subset_taxa(ps_filtered, Family  %in% keepPhyla) 


top10family_mergedv2<-ps_filtered%>%
  subset_taxa(Family %in% top10_family_taxa)%>%
  psmelt()%>%
  as.data.frame()%>%
  group_by(Site,Condition)%>%
  mutate(Abundance = Abundance / sum(Abundance)*100)%>%
  # mutate(Family= factor(Family,levels = top10_family_taxa))%>%
  ggplot(aes(fill=Family,x=Site,y=Abundance))+
  geom_bar(aes(), stat="identity", position="stack")+geom_col(color=NA)+
  theme_classic()+
  theme(legend.key.size = unit(12,"pt"),
        legend.title =  element_text(size = 11,face = "bold"),
        legend.text = element_text(face = "italic"),
        strip.text.x = element_text(face = "bold", size = 13),
        # strip.text.x = element_blank(),
        # strip.background = element_blank(),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_text(face="bold", colour = "black",size = 8),
        axis.ticks.x= ,
        axis.title.y = element_markdown(face="bold", colour = "black",size = 12),
        legend.position = "none",
        plot.title = element_text(size = 11,face = "bold",hjust = 0.5))+
  facet_wrap(~ Condition, scales = "free_x")+
  labs(x=NULL ,y="Relative Abundance (%)", title = "Family")+
  # guides(x = "none")+
  scale_fill_manual(name="Family",values=ggsci::pal_d3(palette = "category10")(10))

top10family_mergedv2

# ggsave("top10family_mergedv2.png",plot = top10family_mergedv2, height = 3.5, width = 7.5)


# merge plots
### supplementary
merged_general <-plot_grid(general_phylum,top10family,ncol = 1,nrow = 2,
                           labels = c("A.1","A.2"))

# merged
merged_merged_plots <-plot_grid(general_phylum_merged,top10family_merged,ncol = 2,nrow = 1,labels = c("B.1","B.2"))


#merged v2 

merged_mergedV2_plots <-plot_grid(general_phylum_mergedv2,top10family_mergedv2,ncol = 2,nrow = 1,labels = c("C.1","C.2"))


final_merged<- plot_grid(merged_general,plot_grid(merged_merged_plots,merged_mergedV2_plots,rel_widths = c(1,1.5)),ncol = 1,nrow = 2,
                         rel_heights = c(1.5,1))
final_merged

ggsave("supplementary_merged_plot_final.tiff",plot = final_merged, height = 10, width = 12)


##### box plots and rel ab comparisons ######

####### General review box plots

# phylum

# testing rel abundance 

abundance_counts<- ps_rm_taxa%>%
  psmelt()%>%
  as.data.frame()

library(nortest)
ad.test(abundance_counts$Abundance)
# p<0.05 so our data is not normally distributed > wilcox | kruskal

phylum_rel_test <- ps_rm_taxa  %>% 
  psmelt()%>%
  group_by(Phylum) %>% 
  wilcox_test(data = ., Abundance ~ Condition) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")%>%
  add_xy_position(scales = "free")

ps_filtered_phylum <- tax_glom(ps_rm_taxa, "Phylum")

phylum_review <- ps_filtered_phylum%>%
  psmelt%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  theme_bw()+
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ Phylum, scales = "free",ncol = 4,nrow = 2)+
  theme(legend.position = "none")+
  stat_pvalue_manual(phylum_rel_test, label = "p.adj.signif")

phylum_review

ggsave("phylum_review_final.tiff",plot = phylum_review, height = 8, width = 10)

# family
# testing rel abundance 
family_rel_test <- ps_rm_taxa  %>% 
  psmelt()%>%
  group_by(Family) %>% 
  wilcox_test(data = ., Abundance ~ Condition) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")%>%
  add_xy_position(scales = "free")%>%
  filter(Family %in% top10_family_taxa)


ps_filtered_family_boxplot <- subset_taxa(ps_rm_taxa, Family %in% top10_family_taxa)

family_review <- psmelt(ps_filtered_family_boxplot)%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Family), height = 0, width = .2) +
  theme_bw()+
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ Family, scales = "free",nrow = 2,ncol = 5)+
  theme(legend.position = "none")+
  stat_pvalue_manual(family_rel_test, label = "p.adj.signif")

family_review

ggsave("family_review_final.tiff",plot = family_review, height = 8, width = 10)


# genus
# testing rel abundance 
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

ps_filtered_genus_boxplot <- subset_taxa(ps_rm_taxa, Genus %in% top10_genus_taxa)

Genus_review <- psmelt(ps_filtered_genus_boxplot)%>%
  ggplot(aes(x = Condition, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  theme_bw()+
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ Genus, scales = "free", nrow = 3,ncol = 5)+
  theme(legend.position = "none")+
  stat_pvalue_manual(genus_rel_test, label = "p.adj.signif")

Genus_review

ggsave("Genus_review_final.tiff",plot = Genus_review, height = 8, width = 10)


###### Alpha Div ##### #####
set.seed(16022000)
alpha_diversity_metrics <- estimate_richness(ps_rarefied, measures = c("Observed", 
                                                                       # "Chao1", 
                                                                       # "ACE",
                                                                       "Shannon", "Simpson"))
metadata <- data.frame(ps_rarefied@sam_data)
metadata <- cbind(metadata, alpha_diversity_metrics) 


# if the results of shapiro test is greater than 0.05, the data is distributed normal (use anova, t.test), if not use (kruskal, wilcox)
shapiro.test(metadata$Observed)
shapiro.test(metadata$Chao1)
shapiro.test(metadata$ACE)
shapiro.test(metadata$Shannon)
shapiro.test(metadata$Simpson)

# Condition
t_test(data = metadata, formula = Observed ~ Condition,p.adjust.method = "fdr") 
t_test(data = metadata, formula = Chao1 ~ Condition,p.adjust.method = "fdr") 
t_test(data = metadata, formula = ACE ~ Condition,p.adjust.method = "fdr") 
wilcox_test(data = metadata, formula = Shannon ~ Condition,p.adjust.method = "fdr")
wilcox_test(data = metadata, formula = Simpson ~ Condition,p.adjust.method = "fdr")

# Site
anova_test(data = metadata, formula = Observed ~ Site) 
anova_test(data = metadata, formula = Chao1 ~ Site) 
anova_test(data = metadata, formula = ACE ~ Site) 
kruskal_test(data = metadata, formula = Shannon ~ Site)
kruskal_test(data = metadata, formula = Simpson ~ Site)

alpha_indexes_plot <-plot_richness(ps_rarefied,x = "Condition", color = "Condition",measures = c("Observed", 
                                                                                                 # "Chao1","ACE", 
                                                                                                 "Shannon", "Simpson"))+
  geom_boxplot()+
  theme_bw()+
  scale_color_lancet(name="Diagnosis")+
  labs(x="Diagnosis")

ggsave(plot = alpha_indexes_plot, "alpha_indexes_plot_final.tiff",height = 4,width = 5.5)