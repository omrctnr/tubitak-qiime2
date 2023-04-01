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
permanova <- adonis2(bc_distance~ Condition, data=metadata,permutations = 9999)
permanova
p.adjust(permanova$`Pr(>F)`, method = "fdr")

## p= 0.043
## R= 0.025

## PCoA Plot 

pcoa_plot <- metadata %>%  
  mutate(Condition = replace(Condition, Condition == "AL", "Aseptic Loosening"),
         Condition = replace(Condition, Condition == "PJI", "Prosthetic Joint Infection")) %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color= Condition), fill= "white", show.legend=FALSE)+
  geom_point() + stat_ellipse(linetype=2) + 
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold",size = 9, colour = "black"),
        axis.text.y = element_text(face = "bold",size = 9, colour = "black"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(20,"pt"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.text = element_text(face="bold"),
        plot.caption = element_markdown(face="bold", size = 11),
        panel.border = element_rect(colour = "black", fill = NA,size=1))+
  scale_color_manual(breaks = c("Aseptic Loosening","Prosthetic Joint Infection"),values=c("blue", "red"))+
  labs(x = "PC1 [8.4%]", y= "PC2 [7.6%]",
       # caption = "<i>p</i><0.05, PERMANOVA"
  )

pcoa_plot

ggsave("pcoa_plot_gg.png",plot=pcoa_plot, height=4.5, width=4.5)  


## ## ## ## ##  
## 3D PCoA  ##
## ## ## ## ## 
library(rgl)
metadata_3d<-metadata %>%
  mutate(condition_color= case_when(Condition == "AL" ~ "blue",
                                    Condition == "PJI" ~ "red",
                                    TRUE ~ NA_character_))

plot3d(x=metadata_3d$pcoa1, y=metadata_3d$pcoa2, z= metadata_3d$pcoa3,
       xlab = "PC1", ylab = "PC2", zlab = "PC3",
       col= metadata_3d$condition_color, type = "s", size = 1)
snapshot3d("gif/3d_pcoa_rgl.png")

## the code below can be used to create an animated gif file 
movie3d(spin3d(axis=c(0,0,1), rpm = 2), duration=10, dir="./")