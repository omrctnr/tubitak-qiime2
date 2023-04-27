library(forcats) # fct_reorder() function
library(microbiomeMarker) # conduct LEfSe
library(dplyr)
library(tidyr)   # filter and reformat data frames
library(stringr) # str_detect() function
set.seed(123123)


## lefse analysis on family, genus and species

ps_rarefied <- rarefy_even_depth(ps_filtered,rngseed=123123, 
                                 sample.size=min(sample_sums(ps_filtered)), replace = FALSE)


lefse_ps_rarefied <-run_lefse(ps_rarefied,
                          wilcoxon_cutoff = 0.05, norm = "none",multigrp_strat = FALSE, kw_cutoff = 0.05,lda_cutoff = 1,
                          group = 'Condition')

marker_table_rarefied <- as_data_frame(lefse_ps_rarefied@marker_table)


## rarefied ps object

lefse_plot_rarefied <- marker_table_rarefied %>%
  
  separate(feature, sep = "\\|", remove = FALSE, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"))%>%
  mutate(feature_main = case_when(str_detect(feature, "s__") ~ str_replace_all(str_extract(feature, "s__.*"), "__", "_"),
                                  !str_detect(feature, "s__") & str_detect(feature, "g__")~ str_replace_all(str_extract(feature, "g__.*"), "__", "_"),
                                  !str_detect(feature, "g__") & str_detect(feature, "f__")~ str_replace_all(str_extract(feature, "f__.*"), "__", "_"),
                                  !str_detect(feature, "f__") & str_detect(feature, "o__")~ str_replace_all(str_extract(feature, "o__.*"), "__", "_"),
                                  !str_detect(feature, "o__") & str_detect(feature, "c__")~ str_replace_all(str_extract(feature, "c__.*"), "__", "_"),
                                  !str_detect(feature, "c__") & str_detect(feature, "p__")~ str_replace_all(str_extract(feature, "p__.*"), "__", "_"),
                                  !str_detect(feature, "p__") & str_detect(feature, "k__")~ str_replace_all(str_extract(feature, "k__.*"), "__", "_"),
                                  TRUE ~ NA_character_),
         enrich_group = replace(enrich_group, enrich_group == "AL", "Aseptic Loosening"),
         enrich_group = replace(enrich_group, enrich_group == "PJI", "Prosthetic Joint Infection"),
         feature_main = replace(feature_main, feature_main == "s_dispar", 's_Veillonella <span class="tab"></span> dispar'),
         feature_main = replace(feature_main, feature_main == "s_parainfluenzae", 's_Haemophilus <span class="tab"></span> parainfluenzae'))%>%
  filter(!feature_main=="s_Streptococcus_s_")%>%
  mutate(feature_main = sub("(.{2})(.*)", "\\1*\\2", feature_main),
         feature_main = paste(feature_main,"*", sep = "")) %>%
  mutate(ef_lda = if_else(enrich_group=="Aseptic Loosening", -1 * ef_lda, ef_lda),
         feature_main = fct_reorder(feature_main, ef_lda))%>%

  ggplot(aes(x=feature_main,y=ef_lda))+
  geom_bar(aes(fill=enrich_group),stat = 'identity', col="black",width = 0.6)+
  labs(y="LDA Score (log 10)", x=NULL, title = NULL)+
  theme_minimal() +
  coord_flip()+
  scale_y_continuous(expression(log[10](italic("LDA Score"))),
                     breaks = seq(-3,3, by=1), limits = c(-2, 2))+
  scale_fill_manual(breaks = c("Aseptic Loosening","Prosthetic Joint Infection"), 
                    values = c("forestgreen", "goldenrod"))+
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
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank(),
        plot.background = element_rect(),
        plot.caption = element_text(face = "italic", size = 12),
  )

lefse_plot_rarefied

################# lefse analysis on site parameter
ps_rarefied <- rarefy_even_depth(ps_filtered,rngseed=123123, 
                                 sample.size=min(sample_sums(ps_filtered)), replace = FALSE)

lefse_ps_rarefied_site <-run_lefse(ps_rarefied,
                                   wilcoxon_cutoff = 0.05, norm = "none",multigrp_strat = FALSE, kw_cutoff = 0.05,lda_cutoff = 1,
                                   group = "Site")

marker_table_rarefied_site <- as_data_frame(lefse_ps_rarefied_site@marker_table)


lefse_plot_rarefied_site<-marker_table_rarefied_site %>%
  
  separate(feature, sep = "\\|", remove = FALSE, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"))%>%
  mutate(feature_main = case_when(str_detect(feature, "s__") ~ str_replace_all(str_extract(feature, "s__.*"), "__", "_"),
                                  !str_detect(feature, "s__") & str_detect(feature, "g__")~ str_replace_all(str_extract(feature, "g__.*"), "__", "_"),
                                  !str_detect(feature, "g__") & str_detect(feature, "f__")~ str_replace_all(str_extract(feature, "f__.*"), "__", "_"),
                                  !str_detect(feature, "f__") & str_detect(feature, "o__")~ str_replace_all(str_extract(feature, "o__.*"), "__", "_"),
                                  !str_detect(feature, "o__") & str_detect(feature, "c__")~ str_replace_all(str_extract(feature, "c__.*"), "__", "_"),
                                  !str_detect(feature, "c__") & str_detect(feature, "p__")~ str_replace_all(str_extract(feature, "p__.*"), "__", "_"),
                                  !str_detect(feature, "p__") & str_detect(feature, "k__")~ str_replace_all(str_extract(feature, "k__.*"), "__", "_"),
                                  TRUE ~ NA_character_),
         enrich_group = replace(enrich_group, enrich_group == "Synovial_Fluid", "Synovial Fluid"))%>%
  filter(!feature_main=="s_Finegoldia_s_")%>%
  filter(!feature_main== "s_Fusobacterium_s_")%>%
  mutate(feature_main = sub("(.{2})(.*)", "\\1*\\2", feature_main),
         feature_main = paste(feature_main,"*", sep = "")) %>%
  mutate(ef_lda = if_else(enrich_group=="Synovial Fluid", -1 * ef_lda, ef_lda),
         feature_main = fct_reorder(feature_main, ef_lda))%>%
  
  ggplot(aes(x=feature_main,y=ef_lda))+
  geom_bar(aes(fill=enrich_group),stat = 'identity', col="black",width = 0.6)+
  labs(y="LDA Score (log 10)", x=NULL, title = NULL)+
  theme_minimal() +
  coord_flip()+
  scale_y_continuous(expression(log[10](italic("LDA Score"))),
                     breaks = seq(-3,3, by=1), limits = c(-2, 2))+
  scale_fill_manual(breaks = c("Synovial Fluid"), 
                    values = c("red"))+
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
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank(),
        plot.background = element_rect(),
        plot.caption = element_text(face = "italic", size = 12),
  )

lefse_plot_rarefied_site

### combine plots
library(cowplot) # merging function
lefse_plot_rarefied
lefse_plot_rarefied_site

d1<- plot_grid(lefse_plot_rarefied,lefse_plot_rarefied_site, nrow = 1,ncol = 2,labels = c("A","B"),
               label_size = 15, rel_widths = c(1.2,1))

d1

ggsave("lefse_plots.png",plot = d1,height = 5, width = 10)
