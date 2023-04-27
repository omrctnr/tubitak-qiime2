library(ggpubr)
library(rstatix)
library(ggprism)
library(forcats)
set.seed(123123)

ps_filt_rel_ab <- transform_sample_counts(ps_filtered, function(x)100* x / sum(x))
ps_filt_rel_ab_melt<-psmelt(ps_filt_rel_ab)

ps_filt_rel_ab_melt<-ps_filt_rel_ab_melt %>%
  mutate(Species = replace(Species, Species == "fragilis", "Bacteroides fragilis"))

# relative abundance testing for species
signif_spec <- ps_filt_rel_ab_melt%>%
  group_by(Species)%>%
  select(Species, Abundance, Site)%>%
  drop_na()%>%
  anova_test(Abundance~Site)%>%
  adjust_pvalue(method = "fdr")%>%
  add_significance("p.adj")%>%
  as_data_frame()%>%
  filter(p.adj < 0.05)

tukey_spec <- ps_filt_rel_ab_melt %>% 
  group_by(Species)%>%
  filter(Species == signif_spec[["Species"]])%>%
  tukey_hsd(Abundance ~ Site) %>% 
  adjust_pvalue(method = "fdr")%>%
  add_significance("p.adj") %>% 
  add_y_position()

plot_signif_species<- ps_filt_rel_ab_melt%>%
  filter(Species == signif_spec[["Species"]])%>%
ggboxplot(x="Site", y="Abundance",palette = "jco", facet.by = "Species")+
  stat_pvalue_manual(tukey_spec, hide.ns = T)


# relative abundance testing for genera

signif_gen <- ps_filt_rel_ab_melt%>%
  group_by(Genus)%>%
  select(Genus, Abundance, Site)%>%
  drop_na()%>%
  anova_test(Abundance~Site)%>%
  adjust_pvalue(method = "fdr")%>%
  add_significance("p.adj")%>%
  as_data_frame()%>%
  filter(p.adj < 0.05)

tukey_gen <- ps_filt_rel_ab_melt %>% 
  group_by(Genus)%>%
  filter(Genus == signif_gen[["Genus"]])%>%
  tukey_hsd(Abundance ~ Site) %>% 
  adjust_pvalue(method = "fdr")%>%
  add_significance("p.adj") %>%
  add_y_position()

plot_signif_genus_total<- ps_filt_rel_ab_melt%>%
  filter(Genus == signif_gen[["Genus"]])%>%
  ggboxplot(x="Site", y="Abundance",palette = "jco", facet.by = "Genus")+
  stat_pvalue_manual(tukey_gen, hide.ns = T)


tukey_gen_1 <- ps_filt_rel_ab_melt %>% 
  group_by(Genus)%>%
  filter(Genus == signif_gen[["Genus"]][1])%>%
  tukey_hsd(Abundance ~ Site) %>% 
  adjust_pvalue(method = "fdr")%>%
  add_significance("p.adj") %>%
  add_y_position()

plot_signif_genus_1<- ps_filt_rel_ab_melt%>%
  filter(Genus == signif_gen[["Genus"]][1])%>%
  ggboxplot(x="Site", y="Abundance",palette = "jco", facet.by = "Genus")+
  stat_pvalue_manual(tukey_gen_1, hide.ns = T)

tukey_gen_2 <- ps_filt_rel_ab_melt %>% 
  group_by(Genus)%>%
  filter(Genus == signif_gen[["Genus"]][2])%>%
  tukey_hsd(Abundance ~ Site) %>% 
  adjust_pvalue(method = "fdr")%>%
  add_significance("p.adj") %>%
  add_y_position()

plot_signif_genus_2<- ps_filt_rel_ab_melt%>%
  filter(Genus == signif_gen[["Genus"]][2])%>%
  mutate(Site = factor(Site, levels = c("Synovial_Fluid", "Biofilm","Tissue")))%>%
  ggboxplot(x="Site", y="Abundance",palette = "jco", facet.by = "Genus")+
  stat_pvalue_manual(tukey_gen_2, hide.ns = T)


# relative abundance testing for families

signif_family <- ps_filt_rel_ab_melt%>%
  group_by(Family)%>%
  select(Family, Abundance, Site)%>%
  drop_na()%>%
  anova_test(Abundance~Site)%>%
  adjust_pvalue(method = "fdr")%>%
  add_significance("p.adj")%>%
  as_data_frame()%>%
  filter(p.adj < 0.05)

tukey_family <- ps_filt_rel_ab_melt %>% 
  group_by(Family)%>%
  filter(Family == signif_family[["Family"]])%>%
  tukey_hsd(Abundance ~ Site) %>% 
  adjust_pvalue(method = "fdr")%>%
  add_significance("p.adj") %>% 
  add_y_position()

plot_signif_family<- ps_filt_rel_ab_melt%>%
  filter(Family == signif_family[["Family"]])%>%
  ggboxplot(x="Site", y="Abundance",palette = "jco", facet.by = "Family")+
  stat_pvalue_manual(tukey_family, hide.ns = T)


library(cowplot) # merging function
plot_signif_species<- plot_signif_species +theme_bw() + labs(x="", y="Relative Abundance (%)") + theme(axis.title.y = element_text(size = 9))
plot_signif_genus_1<- plot_signif_genus_1 +theme_bw() + labs(x="", y="") + theme(axis.title.y = element_text(size = 9))
plot_signif_genus_2<- plot_signif_genus_2 +theme_bw() + labs(x="", y="Relative Abundance (%)") + theme(axis.title.y = element_text(size = 9))
plot_signif_family<- plot_signif_family +theme_bw() + labs(x="", y="") + theme(axis.title.y = element_text(size = 9))

d1<- plot_grid(plot_signif_species,
               plot_signif_genus_1,
               plot_signif_genus_2,
               plot_signif_family, nrow = 2,ncol = 2)

ggsave("sig_taxa_site.png", plot = d1, width = 5,height = 4.1)

# compare_means(Abundance ~ Site,
#                   method = "anova",
#                   group.by = "Species",
#                   psmelt(ps_filt_rel_ab),
#                   p.adjust.method = "fdr")