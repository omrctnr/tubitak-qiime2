sum(subset_samples(ps_filtered, Condition =="Prosthetic Joint Infection")@otu_table)
# 2731506
sum(subset_samples(ps_filtered, Condition =="Aseptic Loosening")@otu_table)
# 50788

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
sort(taxa_sums(ps_sort_Genus),decreasing = TRUE)
sort(taxa_sums(ps_sort_Genus),decreasing = TRUE)[1] / sum(taxa_sums(ps_sort_Genus)) * 100
sort(taxa_sums(ps_sort_Genus),decreasing = TRUE)[2] / sum(taxa_sums(ps_sort_Genus)) * 100
sort(taxa_sums(ps_sort_Genus),decreasing = TRUE)[3] / sum(taxa_sums(ps_sort_Genus)) * 100
sort(taxa_sums(ps_sort_Genus),decreasing = TRUE)[4] / sum(taxa_sums(ps_sort_Genus)) * 100
sort(taxa_sums(ps_sort_Genus),decreasing = TRUE)[5] / sum(taxa_sums(ps_sort_Genus)) * 100

#Species
# ps objesinde varolan specieslerin gruplara göre ayrimi

ps_al<- subset_samples(ps_filtered, Condition == "Aseptic Loosening")
ps_pji<- subset_samples(ps_filtered, Condition == "Prosthetic Joint Infection")

# gruplarda var olan spesifik türlerin bir araya getirilip isimlerinin ASV'ler ile degistirilmesi 
ps_sort_Species_al <- tax_glom(ps_al, "Species")
taxa_names(ps_sort_Species_al) <- tax_table(ps_sort_Species_al)[, 7]

ps_sort_Species_pji <- tax_glom(ps_pji, "Species")
taxa_names(ps_sort_Species_pji) <- tax_table(ps_sort_Species_pji)[, 7]


# iki gruba ait specieslerin ve toplam countslarinin oldugu tablo
al_species_count <- as.data.frame((sort(taxa_sums(ps_sort_Species_al),decreasing = TRUE)))%>%
  rownames_to_column("Species")
colnames(al_species_count)[2] <- "al_read_count"

pji_species_count<- as.data.frame((sort(taxa_sums(ps_sort_Species_pji),decreasing = TRUE)))%>%
  rownames_to_column("Species")
colnames(pji_species_count)[2] <- "pji_read_count"

all_species_count <- full_join(al_species_count, pji_species_count, by = "Species")

# green geneste specieslerin önünde cins ismi olmadigi için
# cins isimlerini belirleyip birlestirmek gerekecek

ps_specific_species <- subset_taxa(ps_filtered, Species %in% taxa_names(ps_sort_Species_al))

specific_species_with_genus <- as.data.frame(tax_table(ps_specific_species)[,6:7]) %>%
  distinct()


final_species <- full_join(all_species_count, specific_species_with_genus, by = "Species") %>%
  unite("full_taxa", c("Genus", "Species"),sep = " ")

# total sayilari ve yüzdeleri hesapla
final_species <- final_species %>%
  arrange(desc(pji_read_count))%>%
  mutate(total = al_read_count + pji_read_count,
         al_percentage = al_read_count/ total * 100,
         pji_percentage = pji_read_count/ total* 100)%>%
  mutate(al_read_count_perc = paste( " (% ",substr(al_percentage,1,4), ")",sep = ""),
         pji_read_count_perc =  paste(" (% ",substr(pji_percentage,1,4), ")",sep = ""))%>%
  select(-al_percentage, -pji_percentage)

# excelde düzenlenmesi için extract edildi
writexl::write_xlsx(final_species, "final_species_table1.xlsx")

# most species by diagnosis group

# pji group

sort(taxa_sums(ps_sort_Species_pji),decreasing = TRUE)
sort(taxa_sums(ps_sort_Species_pji),decreasing = TRUE)[1] / sum(taxa_sums(ps_sort_Species_pji)) * 100
sort(taxa_sums(ps_sort_Species_pji),decreasing = TRUE)[2] / sum(taxa_sums(ps_sort_Species_pji)) * 100
sort(taxa_sums(ps_sort_Species_pji),decreasing = TRUE)[3] / sum(taxa_sums(ps_sort_Species_pji)) * 100
sort(taxa_sums(ps_sort_Species_pji),decreasing = TRUE)[4] / sum(taxa_sums(ps_sort_Species_pji)) * 100
sort(taxa_sums(ps_sort_Species_pji),decreasing = TRUE)[5] / sum(taxa_sums(ps_sort_Species_pji)) * 100

# al group
sort(taxa_sums(ps_sort_Species_al),decreasing = TRUE)
sort(taxa_sums(ps_sort_Species_al),decreasing = TRUE)[1] / sum(taxa_sums(ps_sort_Species_al)) * 100
sort(taxa_sums(ps_sort_Species_al),decreasing = TRUE)[2] / sum(taxa_sums(ps_sort_Species_al)) * 100
sort(taxa_sums(ps_sort_Species_al),decreasing = TRUE)[3] / sum(taxa_sums(ps_sort_Species_al)) * 100
sort(taxa_sums(ps_sort_Species_al),decreasing = TRUE)[4] / sum(taxa_sums(ps_sort_Species_al)) * 100
sort(taxa_sums(ps_sort_Species_al),decreasing = TRUE)[5] / sum(taxa_sums(ps_sort_Species_al)) * 100


# library(kableExtra)
# 
# final_species %>%
#   arrange(desc(pji_read_count))%>%
#   mutate(total = al_read_count + pji_read_count,
#          al_percentage = al_read_count/ total * 100,
#          pji_percentage = pji_read_count/ total* 100)%>%
#   # mutate_at(c('al_percentage', 'pji_percentage'), as.numeric)%>%
#   mutate(al_read_count_perc = paste( " (% ",substr(al_percentage,1,4), ")",sep = ""),
#          pji_read_count_perc =  paste(" (% ",substr(pji_percentage,1,4), ")",sep = ""))%>%
#   select(-al_percentage, -pji_percentage) %>%
#   kbl() %>%
#   kable_paper()
# # add_header_above(c(" " = 1, "Read Counts" = 3, " "= 2))

