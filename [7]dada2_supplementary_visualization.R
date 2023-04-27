## ## ## ## ## ## ##
## VISUALIZATION  ##
## ## ## ## ## ## ##

## general review phyla with stacked bar plot strafied by condition and site

ps_filt_rel_ab <- transform_sample_counts(ps_filtered, function(x)100* x / sum(x))

general_phylum_site <- ps_filt_rel_ab %>%
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
  facet_wrap(Condition~ Site, scales = "free_x",nrow = 1,ncol = 5)+
  labs(name="Phylum",x=NULL ,y="   ", title = NULL)+
  guides(x = "none")+
  # scale_fill_manual(name=NULL,values=c(brewer.pal(9,"Paired"),"gray"))
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#33A02C",
                               "#B2DF8A","#E31A1C","#FB9A99",
                               "#FDBF6F","#FF7F00","#CAB2D6"))

general_phylum_site

## top 10 family with stacked bar plot strafied by Condition

# select top 10 ASVs in a phyloseq object
top10_family <- top_taxa(ps_filtered, tax_level = "Family", 
                         n_taxa = 10,include_na_taxa = FALSE)

df_top10_family <- top10_family$top_taxa
vec_top10_family <- df_top10_family[["Family"]]

ps_filtered_family <- subset_taxa(ps_filtered, Family %in% vec_top10_family)

ps_filtered_family <- transform_sample_counts(ps_filtered_family, function(x)100* x/sum(x))

top10family_site<- plot_bar(ps_filtered_family, fill= "Family")+
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
  facet_wrap(Condition~Site, scales = "free_x",nrow = 1,ncol = 5)+
  labs(x=NULL ,y="Relative Abundance <br>(%)", title = NULL)+
  guides(x = "none")+
  scale_fill_manual(name="Family",values=c(brewer.pal(10,"Paired"),"gray"))

top10family_site

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

top15Genus_site<- plot_bar(ps_filtered_Genus, fill= "Genus")+
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
  facet_wrap(Condition~Site, scales = "free_x",nrow = 1,ncol = 5)+
  labs(x=NULL ,y="   ", title = NULL)+
  guides(x = "none")+
  scale_fill_manual(name="Genus",values = rev(distinct_palette(n = 16, pal = "kelly",add = NA)))

top15Genus_site

ggsave("top15Genus_gg.png",plot = top15Genus, height = 3.5, width = 8.5)

##### merge plots
library(cowplot) # merging function

general_phylum_site
top10family_site
top15Genus_site

d1<- plot_grid(general_phylum_site,top10family_site,top15Genus_site, nrow = 3,ncol = 1,labels = c("A","B","C"),
               label_size = 15)


ggsave("merged_site_test.tiff", plot = d1, height=8.5, width=15)

