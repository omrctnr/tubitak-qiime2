library(readxl)
library(phyloseq)
library(ggplot2)      # graphics
library(RColorBrewer) # Special Color Palettes
library(fantaxtic)    # get top n taxa
library(dplyr)  # dataframe manipulation
library(ggtext) # allows writing HTML code to change the style of text 
                    # and adding significance labels and stars


mag_mat<- read_excel("qiime_R.xlsx", sheet = "reads")
tax_mat<- read_excel("qiime_R.xlsx", sheet = "TAX")
samples_df <- read_excel("qiime_R.xlsx", sheet = "MET")

mag_mat <- mag_mat %>%
  tibble::column_to_rownames("taxon") 

tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("taxon")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("PatientID") 

mag_mat <- as.matrix(mag_mat)
tax_mat <- as.matrix(tax_mat)

MAG = otu_table(mag_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

carbom <- phyloseq(MAG, TAX, samples)
carbom


## ## ## ## ## ## ## ##
## Stacked Bar Chart ##
## ## ## ## ## ## ## ##

sample_data(carbom)$dummy_var <- -1

# ONLY CONDITION
carbom_top150<-get_top_taxa(carbom,150)
carbom_top150 <- merge_samples(carbom_top150, "Condition")
carbom_top150 = transform_sample_counts(carbom_top150, function(x)100* x/sum(x))

# plot
condition_bar_chart<- plot_bar(carbom_top150, fill= "Phylum")+
  geom_bar(aes( colour = NULL), stat="identity", position="stack")+
  theme_classic()+
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(14,"pt"),
        legend.title = element_text(face = "bold"),
        title=element_text(face = "bold",size = 13),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_text(face="bold", colour = "black",size = 12),
        axis.title.y = element_text(face="bold", colour = "black",size = 12),
        plot.caption.position = "plot",
        plot.caption = element_text(face = "italic",size = 10))+
  labs(x=NULL,y="Relaltive Abundance (%)",
       title = "Phylum Taxonomy Distribution \nAccording to the Prothesis Condition",
       caption = "Created using top 150 OTUs")+
  scale_fill_manual(name=NULL,
                 values=paletteer_d("ggthemes::Classic_Cyclic"))

ggsave("condition_bar_chart.tiff",plot = condition_bar_chart, height = 5, width = 5)

  
# sadece  GENEL BAKIS
carbom_relative = transform_sample_counts(carbom, function(x)100* x/sum(x))

general_phylum <- plot_bar(carbom_relative, fill= "Phylum")+
  geom_bar(aes( colour = Phylum, fill = Phylum ), stat="identity", position="stack")+
  theme_classic()+
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(12,"pt"),
        legend.title = element_blank(),
        title=element_text(face = "bold",size = 13),
        axis.text.y = element_text(face="bold", colour = "black",size = 10),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", colour = "black",size = 12))+
  labs(x="Samples",y="Relaltive Abundance (%)",
       title = "Phylum Taxonomy Distribution")

ggsave("general_phylum.tiff",plot = general_phylum, height = 5, width = 6)

## ## ## ##
## LEfsE ##
## ## ## ##
library(forcats) # fct_reorder() function
library(microbiomeMarker) # calculate LEfSe

lefse_carbom <-run_lefse(carbom,
                  wilcoxon_cutoff = 0.05, norm = "TMM",multigrp_strat = FALSE,kw_cutoff = 0.05,lda_cutoff = 2,
                  group = 'Condition')

marker_table <- as_data_frame(lefse_carbom@marker_table)

marker_table@.Data[[1]][1] <- "f__Bacteroidaceae"
marker_table@.Data[[1]][2] <- "g__Bacteroides"

lefse_plot <- marker_table%>%
                      mutate(ef_lda = if_else(enrich_group=="Infection", -1 * ef_lda, ef_lda),
                             feature = fct_reorder(feature, ef_lda))%>%
                      ggplot(aes(x=feature,y=ef_lda))+
                      geom_bar(aes(fill=enrich_group),stat = 'identity', col="black",width = 0.6)+
                      labs(y="LDA Score (log 10)", x=NULL, title = "LEfSe Analysis")+
                      theme_minimal() +
                      coord_flip()+
                      scale_y_continuous(expression(log[10](italic("LDA score"))),
                                         breaks = seq(-6,6, by=2), limits = c(-6, 6))+
                      scale_fill_manual(breaks = "Infection", values = "red")+
                      geom_text(aes(y = 0, label = feature, hjust = ifelse(ef_lda < 0, -0.03, 1.03)),
                                size=4)+
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

ggsave("lefse_plot.tiff",plot = lefse_plot, height = 5, width = 5)





## ## ## ## ## ## ## ##
##  Alpha Diversity  ##
## ## ## ## ## ## ## ##
library(tidyr)  # dataframe manipulation
library(tibble) # dataframe manipulation


samples_df <- bind_cols(samples_df,estimate_richness(carbom, measures = c("Observed", "Chao1", "ACE", 
                                       "Shannon", "Simpson", "InvSimpson")))
          
wilcox.test(Observed ~ Condition, data=samples_df)
wilcox.test(Chao1 ~ Condition, data=samples_df)
wilcox.test(ACE ~ Condition, data=samples_df)
wilcox.test(se.ACE ~ Condition, data=samples_df)
wilcox.test(Simpson ~ Condition, data=samples_df)
wilcox.test(InvSimpson ~ Condition, data=samples_df)
wilcox.test(Shannon ~ Condition, data=samples_df)


alpha_div_plot <- samples_df %>%
                    mutate(SampleID = rownames(deneme)) %>%
                    gather(key   = alphadiv_index,
                           value = obs_values,
                           -SampleID, -Condition,
                           -ACE, -se.ACE, -Simpson,
                           -Chao1, -se.chao1)%>%
                    mutate(alphadiv_index = factor(alphadiv_index, #sort the index
                                                            levels = c("Observed", "Shannon","InvSimpson")))%>%
                  ggplot(aes(x=Condition, y=obs_values,fill=Condition))+
                    geom_boxplot(show.legend=TRUE, outlier.shape = NA)+
                    theme(strip.background = element_rect(color = NA),
                          strip.text = element_markdown(size = 12, face = "bold", colour = "black"  ),
                          axis.text.x = element_blank(),
                          axis.text.y = element_text(face = "bold",size = 9, colour = "black"),
                          plot.title = element_markdown(hjust = 0.5 ,face="bold", size = 15),
                          axis.title.y = element_text(face = "bold"),
                          legend.title = element_blank(),
                          legend.direction = "horizontal",
                          legend.position = "top",
                          legend.key.size = unit(20,"pt"),
                          legend.key = element_rect(fill = NA, color = NA),
                          legend.text = element_text(face="bold"),
                          plot.caption = element_markdown(face="bold", size = 11))+
                    facet_wrap(~alphadiv_index, scales = "free_y",strip.position = "top")+
                    scale_y_continuous(limits = c(0,NA))+
                    scale_fill_manual(breaks = c("Aseptic Loosing","Infection"),values=c("red", "forestgreen"))+
                    labs(x=NULL, y="Alpha Diversity Measure", caption = "*<i>p</i><0.05, Wilcoxon Test <br><i>ns</i>, no significance")


## add the significance lines and stars

lines<- tibble(
  alphadiv_index = c("Observed", "Shannon", "InvSimpson"),
  
  x = c(0.9,0.9,0.9),
  xend= c(2.1,2.1,2.1),
  y= c(41,3.52,20.5),
  yend= y
)

stars<- tibble(
  alphadiv_index = c("Observed", "Shannon", "InvSimpson"),
  
  x = c(1.5,1.5,1.5),
  y= c(43.5,3.66,21),
  label = c("ns","ns","*")
)

alpha_div_plot <- alpha_div_plot + 
            geom_segment(data= lines, aes(x=x,xend=xend, y=y, yend=yend), 
                         inherit.aes = FALSE,
                         size=0.8)+
            geom_text(data=stars, aes(x=x, y=y, label=label), size=6, inherit.aes = FALSE)

                
ggsave("alpha_div_plot.tiff",plot = alpha_div_plot, height = 5, width = 5)


## ## ## ## ## ## ## 
## Beta Diversity ##
## ## ## ## ## ## ## 
set.seed(123)
library(vegan)  # diversity estimates
library(ape)    #pcoa() function


data_reads <- t(read_excel("qiime_R.xlsx", sheet = "reads"))
colnames(data_reads)<- data_reads[1,]
data_reads<- as.data.frame(data_reads[-1,])

metadata <- read_excel("qiime_R.xlsx", sheet = "MET") %>%
  column_to_rownames(var = "PatientID")

data_reads<-t(dplyr::select(as.data.frame(t(data_reads)),matches(row.names(metadata))))
class(data_reads) <- "numeric"

## create Bray-Curtis Distance
bc_distance <- vegdist(data_reads, method = "bray")
bc_pcoa <- pcoa(bc_distance)

## get principal coordinate values

metadata$pcoa1 <- bc_pcoa$vectors[,1]
metadata$pcoa2 <- bc_pcoa$vectors[,2]
metadata$pcoa3 <- bc_pcoa$vectors[,3]

## get variances of principal coordinate values 
bc_pcoa$values$Relative_eig[[1]] * 100
bc_pcoa$values$Relative_eig[[2]] * 100


## PERMANOVA Test according to the Condition was not significant
adonis2(bc_distance~ Condition, data=metadata,permutations = 999)

## PCoA plot

### alpha grafipinde ldugu gibi y ve x eksenindeki yazilari kalin yaz!

pcoa_plot <-  metadata %>%  
            ggplot(aes(x=pcoa1, y=pcoa2, color= Condition), show.legend=FALSE)+
            geom_point() + stat_ellipse(linetype=2) + 
            theme_classic()+
            theme(plot.title = element_markdown(face="bold", size = 15),
                  axis.text.x = element_text(face = "bold",size = 9, colour = "black"),
                  axis.text.y = element_text(face = "bold",size = 9, colour = "black"),
                  axis.title.y = element_text(face = "bold"),
                  axis.title.x = element_text(face = "bold"),
                  legend.title = element_blank(),
                  legend.position = "top",
                  legend.key.size = unit(20,"pt"),
                  legend.key = element_rect(fill = NA, color = NA),
                  legend.text = element_text(face="bold"),
                  plot.caption = element_markdown(face="bold", size = 11))+
      scale_color_manual(breaks = c("Aseptic Loosing","Infection"),values=c("blue", "red"))+
            labs(x = "PC1 [9.5%]", y= "PC2 [7.1%]",
                 title="Bray-Curtis Distance",
                 caption = "<i>p</i>>0.05, PERMANOVA")

ggsave("pcoa_plot.tiff",plot=pcoa_plot, height=4.5, width=4.5)  

## ## ## ## ##  
## 3D PCoA  ##
## ## ## ## ## 
library(rgl)
metadata_3d<-metadata %>%
  mutate(condition_color= case_when(Condition == "Aseptic Loosing" ~ "blue",
                                   Condition == "Infection" ~ "red",
                                   TRUE ~ NA_character_))

plot3d(x=metadata_3d$pcoa1, y=metadata_3d$pcoa2, z= metadata_3d$pcoa3,
       xlab = "PC1", ylab = "PC2", zlab = "PC3",
       col= metadata_3d$condition_color, type = "s", size = 1)
snapshot3d("3d_pcoa_rgl.png")

## the code below can be used to create an animated gif file 
movie3d(spin3d(axis=c(0,0,1), rpm = 2), duration=10, dir="./")

