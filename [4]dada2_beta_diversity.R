set.seed(123123)
library(vegan)  # diversity estimates
library(ape)    #pcoa() function
library(ggtext)
library(ggplot2)

counts_df <- as.data.frame(ps_filtered@otu_table)

metadata <- read.table("C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/metadata.txt",
                       header = TRUE, sep = "")

## create Bray-Curtis Distance
bc_distance <- vegdist(sqrt(counts_df), method = "bray")
bc_pcoa <- pcoa(bc_distance)

## get principal coordinate values

metadata$pcoa1 <- bc_pcoa$vectors[,1]
metadata$pcoa2 <- bc_pcoa$vectors[,2]
metadata$pcoa3 <- bc_pcoa$vectors[,3]

## get variances of principal coordinate values 
bc_pcoa$values$Relative_eig[[1]] * 100
bc_pcoa$values$Relative_eig[[2]] * 100


## PERMANOVA Test according to the Condition was not significant
permanova_cond <- adonis2(bc_distance~ Condition, data=metadata,permutations = 9999)
permanova_cond
p.adjust(permanova_cond$`Pr(>F)`, method = "fdr")

## p= 0.043
## R= 0.025


permanova_site <- adonis2(bc_distance~ Site, data=metadata,permutations = 9999)
permanova_site

## p= 0.33
## R= 0.039

## PCoA Plot 

pcoa_plot <- metadata %>%  
  mutate(Condition = replace(Condition, Condition == "AL", "Aseptic Loosening"),
         Condition = replace(Condition, Condition == "PJI", "Prosthetic Joint Infection"),
         Site = replace(Site, Site == "Synovial_Fluid", "Synovial Fluid"),
         Site = replace(Site, Site == "Implant", "Biofilm")) %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color= Condition), fill= "white", show.legend=FALSE)+
  geom_point(aes(shape = Site), size = 2) + stat_ellipse(linetype=2) + 
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold",size = 9, colour = "black"),
        axis.text.y = element_text(face = "bold",size = 9, colour = "black"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        # legend.position = "top",
        legend.box = "vertical",
        legend.key.size = unit(20,"pt"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.text = element_text(face="bold"),
        plot.caption = element_markdown(face="bold", size = 11),
        panel.border = element_rect(colour = "black", fill = NA,size=1)
        ) +
  scale_color_manual(breaks = c("Aseptic Loosening","Prosthetic Joint Infection"),values=c("blue", "red"))+
  labs(x = "PC1 [8.4%]", y= "PC2 [7.6%]")+
  guides(color = guide_legend(title = "Diagnosis"),
         shape = guide_legend(title = "Sample Collection Site"))

pcoa_plot

ggsave("pcoa_plot_gg.png",plot=pcoa_plot, height=4.5, width=9.5)  
