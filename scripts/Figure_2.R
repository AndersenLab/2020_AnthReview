library(tidyverse)
library(ggthemes)
library(ggplot2)

##Setwd to folder location

setwd("~/Desktop/2020_AnthReview/")
##Map Figure 1A
world <- map_data("world")
world <- world[world$region != "Antarctica",]

isolation_info <- readr::read_tsv("data/CelegansStrainData.tsv")%>%
  dplyr::filter(isotype_ref_strain == TRUE)%>%
  dplyr::filter(!(release %in% c("20180413", "20200815")))%>%
  dplyr::select(strain, long=longitude, lat=latitude, landscape, substrate)%>%
  dplyr::mutate(n2 = case_when( strain == "N2" ~ "N2",
                                strain != "N2" ~ "None"))


isolation_info$long <- as.numeric(isolation_info$long)
isolation_info$lat <- as.numeric(isolation_info$lat)

map <-ggplot()+ geom_map(data=world, map=world,
                         aes(x=long, y=lat, map_id=region),
                         color="white", fill="#7f7f7f", 
                         size=0.1, alpha=0.5)+
  theme_map()+ 
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

mapfix <-map +geom_point(data = isolation_info%>%filter(n2=="None"),aes(x=long, y=lat, color = n2),size = 0.5)+
  geom_point(data = isolation_info%>%filter(n2=="N2"),aes(x=long, y=lat, color =n2 ),size = 0.5)+
  scale_color_manual(values = c("None" = "black","N2"="Orange"))


##Phenotype change Figure 1B
filtered_data <- rio::import("data/figure1data.tsv")

abz_response <- filtered_data %>%
  dplyr::select(strain,Albendazole_q90.TOF)%>%
  dplyr::filter(!is.na(Albendazole_q90.TOF))
abz_full <- abz_response %>%
  dplyr::mutate( q90_norm = (Albendazole_q90.TOF - min(Albendazole_q90.TOF)) / (max(Albendazole_q90.TOF) - min(Albendazole_q90.TOF)))%>%
  dplyr::mutate(Strain = case_when(strain=="N2" ~ "N2",
                                   TRUE ~ "None"))


response_variation_abz <- abz_full %>%
  ggplot()+
  aes(y =q90_norm)+
  geom_col(aes(x=reorder(strain, Albendazole_q90.TOF),fill=Strain),width = 1)+
  theme_classic(24)+
  ylab("Normalized ABZ Response")+
  scale_y_continuous("Albendazole Response",expand = c(0, 0))+
  xlab("Strain")+
  scale_fill_manual(values = c("N2"="Orange","None"="grey"))+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  theme(axis.title.x=element_text(size=10, vjust=5),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.line.y = element_line(size=0.1),
        axis.line.x = element_line(size=0.1),
        plot.margin = unit(c(0, 0.1, 0, 0.8), "cm"),
        legend.position = "none")

##Phenotype change Figure 1C IVM
ivm_response <- filtered_data %>%
  dplyr::select(strain,Ivermectin_mean.TOF)%>%
  dplyr::filter(!is.na(Ivermectin_mean.TOF))
var_ivm_norm <- ivm_response %>%
  dplyr::mutate( mean_norm = (Ivermectin_mean.TOF - min(Ivermectin_mean.TOF)) / (max(Ivermectin_mean.TOF) - min(Ivermectin_mean.TOF)))

full_ivm <- var_ivm_norm%>%
  dplyr::mutate(Strain = case_when(strain=="N2" ~ "N2",
                                      TRUE ~ "None"))

response_variation_ivm <- full_ivm %>%
  ggplot()+
  aes(y = mean_norm)+
  geom_col(aes(x=reorder(strain,Ivermectin_mean.TOF),fill=Strain),width = 1)+
  theme_classic(24)+
  scale_y_continuous("Ivermectin Response",expand = c(0, 0))+
  xlab("Strain")+
  scale_fill_manual(values = c("N2"="Orange","None"="grey"))+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  theme(axis.title.x=element_text(size=10, vjust=5),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.line.y = element_line(size=0.1),
        axis.line.x = element_line(size=0.1),
        plot.margin = unit(c(0, 0.1, 0, 0.8), "cm"),
        legend.position = "none")

#Figure 1D levamasole 
lev_response <- filtered_data %>%
  dplyr::select(strain,Levamisole_mean.TOF)%>%
  dplyr::filter(!is.na(Levamisole_mean.TOF))
var_lev_norm <- lev_response %>%
  dplyr::mutate( mean_norm = (Levamisole_mean.TOF - min(Levamisole_mean.TOF)) / (max(Levamisole_mean.TOF) - min(Levamisole_mean.TOF)))

full_lev <- var_lev_norm%>%
  dplyr::mutate(Strain = case_when(strain=="N2" ~ "N2",
                                   TRUE ~ "None"))

response_variation_lev <- full_lev %>%
  ggplot()+
  aes(y = mean_norm)+
  geom_col(aes(x=reorder(strain,Levamisole_mean.TOF),fill=Strain),width = 1)+
  theme_classic(24)+
  scale_y_continuous("Levamisole Response",expand = c(0, 0))+
  xlab("Strain")+
  scale_fill_manual(values = c("N2"="Orange","None"="grey"))+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  theme(axis.title.x=element_text(size=10, vjust=5),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.line.y = element_line(size=0.1),
        axis.line.x = element_line(size=0.1),
        plot.margin = unit(c(0, 0.1, 0, 0.8), "cm"),
        legend.position = "none")

figure_2 <- cowplot::plot_grid(mapfix,response_variation_abz,response_variation_ivm,response_variation_lev, nrow = 4, ncol=1,labels = c("A","B","C","D"))
ggsave(filename = "plots/Review_figure2.png", plot = figure_2, device = "png",height = 8,width = 6,units = 'in')
ggsave(filename = "plots/Review_figure2.pdf", plot = figure_2, device = "pdf",height = 8,width = 6,units = 'in')

