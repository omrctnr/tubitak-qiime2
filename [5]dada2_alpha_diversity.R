set.seed(123123)
library(ggtext)
library(ggplot2)
library(phyloseq)

# the rarefaction depth chosen is the minimum sample depth 
# (in this case 604 reads per sample)

ps_rarefied <- rarefy_even_depth(ps_filtered,rngseed=123123, sample.size=min(sample_sums(ps_filtered)), replace = FALSE)


# Calculate diversity for every samples
alpha_diversity_metrics <- estimate_richness(ps_rarefied, measures = c("Observed", "Chao1", 
                                                           "ACE", "Shannon", "Simpson", 
                                                           "InvSimpson"))

metadata <- read.table("C:/Users/ASUS/Desktop/prosthesis_microbiota_study/data/metadata.txt",
                       header = TRUE, sep = "")

metadata <- as.data.frame(sample_data(ps_rarefied))
# Assign the estimated diversity to sample metadata
metadata <- cbind(metadata, alpha_diversity_metrics) 

library(rstatix)

wilcox_test(Observed ~ Condition, data=metadata, p.adjust.method = "fdr")
wilcox_test(Chao1 ~ Condition, data=metadata, p.adjust.method = "fdr")
wilcox_test(ACE ~ Condition, data=metadata, p.adjust.method = "fdr")
wilcox_test(se.ACE ~ Condition, data=metadata, p.adjust.method = "fdr")
wilcox_test(Simpson ~ Condition, data=metadata, p.adjust.method = "fdr")
wilcox_test(InvSimpson ~ Condition, data=metadata, p.adjust.method = "fdr")
wilcox_test(Shannon ~ Condition, data=metadata, p.adjust.method = "fdr")

summary(aov(data = metadata, Observed ~ Site ))
summary(aov(data = metadata, Chao1 ~ Site ))
summary(aov(data = metadata, ACE ~ Site ))
summary(aov(data = metadata, se.ACE ~ Site ))
summary(aov(data = metadata, Simpson ~ Site ))
summary(aov(data = metadata, InvSimpson ~ Site ))
summary(aov(data = metadata, Shannon ~ Site ))

## no significant any value

alpha_div_plot<- metadata %>%
  mutate(SampleID = rownames(metadata)) %>%
  gather(key   = alphadiv_index,
         value = obs_values,
         -SampleID, -Condition,
         -ACE, -se.ACE, -Simpson,
         -Chao1, -se.chao1)%>%
  mutate(alphadiv_index = factor(alphadiv_index, #sort the index
                                 levels = c("Observed", "Shannon","InvSimpson")),
         Condition = replace(Condition, Condition == "AL", "Aseptic Loosening"),
         Condition = replace(Condition, Condition == "PJI", "Prosthetic Joint Infection"))%>%
  ggplot(aes(x=Condition, y=obs_values,fill=Condition))+
  geom_boxplot(show.legend=TRUE, outlier.shape = NA)+
  theme_minimal()+
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
        plot.caption = element_markdown(face="bold", size = 11),
        panel.border = element_rect(colour = "black", fill = NA,size=1))+
  facet_wrap(~alphadiv_index, scales = "free_y",strip.position = "top")+
  scale_y_continuous(limits = c(0,NA))+
  scale_fill_manual(breaks = c("Aseptic Loosening","Prosthetic Joint Infection"),values=c("blue", "red"))+
  labs(x=NULL, y="Alpha Diversity Metrics")

alpha_div_plot

ggsave("alpha_div_plot_gg.png",plot = alpha_div_plot, height = 5, width = 5)
