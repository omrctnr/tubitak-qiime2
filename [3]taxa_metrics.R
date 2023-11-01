
sum(subset_samples(ps_filtered, Condition =="PJI")@otu_table)
# 2866731
sum(subset_samples(ps_filtered, Condition =="AL")@otu_table)
# 64710

2866731/sum(ps_filtered@otu_table)

summary(sample_sums(ps_filtered))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 604    1864    3034   51524    4902  961041 


tax_glom(ps_filtered, "Phylum")
tax_glom(ps_filtered, "Family")
tax_glom(ps_filtered, "Genus")
tax_glom(ps_filtered, "Species")

ps_phylum <- tax_glom(ps_filtered, "Phylum")
taxa_names(ps_phylum) <- tax_table(ps_phylum)[, 2]
sort(taxa_sums(ps_phylum),decreasing = TRUE)
sort(taxa_sums(ps_phylum),decreasing = TRUE)[1] / sum(taxa_sums(ps_phylum)) * 100
sort(taxa_sums(ps_phylum),decreasing = TRUE)[2] / sum(taxa_sums(ps_phylum)) * 100
sort(taxa_sums(ps_phylum),decreasing = TRUE)[3] / sum(taxa_sums(ps_phylum)) * 100
sort(taxa_sums(ps_phylum),decreasing = TRUE)[4] / sum(taxa_sums(ps_phylum)) * 100
sort(taxa_sums(ps_phylum),decreasing = TRUE)[5] / sum(taxa_sums(ps_phylum)) * 100

ps_sort_Family <- tax_glom(ps_filtered, "Family")
taxa_names(ps_sort_Family) <- tax_table(ps_sort_Family)[, 5]
sort(taxa_sums(ps_sort_Family),decreasing = TRUE)
sort(taxa_sums(ps_sort_Family),decreasing = TRUE)[1] / sum(taxa_sums(ps_sort_Family)) * 100
sort(taxa_sums(ps_sort_Family),decreasing = TRUE)[2] / sum(taxa_sums(ps_sort_Family)) * 100
sort(taxa_sums(ps_sort_Family),decreasing = TRUE)[3] / sum(taxa_sums(ps_sort_Family)) * 100
sort(taxa_sums(ps_sort_Family),decreasing = TRUE)[4] / sum(taxa_sums(ps_sort_Family)) * 100
sort(taxa_sums(ps_sort_Family),decreasing = TRUE)[5] / sum(taxa_sums(ps_sort_Family)) * 100

ps_sort_Genus <- tax_glom(ps_filtered, "Genus")
taxa_names(ps_sort_Genus) <- tax_table(ps_sort_Genus)[, 6]

ps_sort_Genus <- tax_glom(subset_samples(ps_filtered, Condition =="AL"), "Genus")
taxa_names(ps_sort_Genus) <- tax_table(ps_sort_Genus)[, 6]
sort(taxa_sums(ps_sort_Genus),decreasing = TRUE) /sum(taxa_sums(ps_sort_Genus)) * 100


#Species
# ps objesinde varolan specieslerin gruplara g?re ayrimi

ps_al<- subset_samples(ps_filtered, Condition == "AL")
ps_pji<- subset_samples(ps_filtered, Condition == "PJI")

# gruplarda var olan spesifik t?rlerin bir araya getirilip isimlerinin ASV'ler ile degistirilmesi 
ps_sort_Species_al <- tax_glom(ps_al, "Species")
taxa_names(ps_sort_Species_al) <- tax_table(ps_sort_Species_al)[, 7]

ps_sort_Species_pji <- tax_glom(ps_pji, "Species")
taxa_names(ps_sort_Species_pji) <- tax_table(ps_sort_Species_pji)[, 7]

asv_species_df<-as.data.frame(tax_table(ps_filtered)[,6:7])

# iki gruba ait specieslerin ve toplam countslarinin oldugu tablo
al_species_count <- as.data.frame((sort(taxa_sums(ps_sort_Species_al),decreasing = TRUE)))%>%
  rownames_to_column("Species")
colnames(al_species_count)[2] <- "al_read_count"

# 1,2,3,4,9

pji_species_count<- as.data.frame((sort(taxa_sums(ps_sort_Species_pji),decreasing = TRUE)))%>%
  rownames_to_column("Species")
colnames(pji_species_count)[2] <- "pji_read_count"

# 1,21,13,20,29

all_species_count <- full_join(al_species_count, pji_species_count, by = "Species")

nelabu<-ps_filtered%>%
  tax_glom("Species")%>%
  subset_samples(Condition=="PJI")%>%
  psmelt()%>%
  as.data.frame()%>%
  group_by(Genus,Species,OTU)%>%
  summarise(Abundance = sum(Abundance))
select(Abundance, Condition, Genus, OTU, Species )


### benzersiz asvleri ve speciesleri tespit et !!
taxa_names_df<-ps_filtered%>%
  tax_table()%>%
  as.data.frame()

taxa_names(ps_filtered) <- paste0(rownames(taxa_names_df)," ",taxa_names_df[["Genus"]]," ",taxa_names_df[["Species"]])
# paste(rownames(taxa_names_df),taxa_names_df[["Genus"]],taxa_names_df[["Species"]])



ps_filt_species<-tax_glom(ps_filtered,"Species")

ps_filt_species_al<- subset_samples(ps_filt_species, Condition =="AL")
ps_filt_species_al_unique<-filter_taxa(ps_filt_species_al, function(x) sum(x) > 0, TRUE)

ps_filt_species_pjı<- subset_samples(ps_filt_species, Condition =="PJI")
ps_filt_species_pjı_unique<-filter_taxa(ps_filt_species_pjı, function(x) sum(x) > 0, TRUE)

all_al_asv<-names(taxa_sums(ps_filt_species_al_unique))
all_pjı_asv<-names(taxa_sums(ps_filt_species_pjı_unique))

names(taxa_sums(ps_filt_species))

unique_for_AL<-setdiff(all_al_asv,all_pjı_asv)
unique_for_PJI<-setdiff(all_pjı_asv,all_al_asv)

al_unique_species_abd<-taxa_sums(ps_filt_species_al_unique)[names(taxa_sums(ps_filt_species_al_unique)) %in% unique_for_AL]
pjı_unique_species_abd<-taxa_sums(ps_filt_species_pjı_unique)[names(taxa_sums(ps_filt_species_pjı_unique)) %in% unique_for_PJI]

cbind(
  as.data.frame(al_unique_species_abd)%>%
  rownames_to_column("Species_AL"),
  as.data.frame(pjı_unique_species_abd)%>%
  rownames_to_column("Species_PJI")
)

excel_workbook<- list("AL_Unique_Species" =   as.data.frame(al_unique_species_abd)%>%
                                                    rownames_to_column("Species_AL"), 
                      "PJI_Unique_Species" =   as.data.frame(pjı_unique_species_abd)%>%
                                                    rownames_to_column("Species_PJI")
                      )

openxlsx::write.xlsx(excel_workbook, "unique_species_by_group.xlsx")




ps_filt_species<-tax_glom(ps_filtered,"Species")

ps_filt_species_al<- subset_samples(ps_filt_species, Condition =="PJI")

sort(taxa_sums(ps_filt_species_al)/ sum(taxa_sums(ps_filt_species_al)) * 100,decreasing = TRUE)

ps_filt_species<-tax_glom(ps_filtered,"Species")

ps_filt_species_al<- subset_samples(ps_filt_species, Condition =="AL")

sort(taxa_sums(ps_filt_species_al)/ sum(taxa_sums(ps_filt_species_al)) * 100,decreasing = TRUE)

