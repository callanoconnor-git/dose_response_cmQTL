library(tidyverse)


nod<- readr::read_csv("~/Downloads/NOD_.75_0_arsenic.csv")
nzo<- readr::read_csv("~/Downloads/NZO_.75_0_arsenic.csv")

dim(nod)

nod<- nod[,c(4,9:15,25:27,30:35)] %>% janitor::clean_names()
nzo<- nzo[,c(4,9:15,25:27,30:35)]%>% janitor::clean_names()
nzo<- nzo %>% dplyr::rename(ar_exposure= "abundance_ratio_log2_0_75_u_m_arsenic_nzo_0_u_m_arsenic_nzo",
                            ar_0 = "abundance_ratio_log2_0_u_m_arsenic_nzo_0_u_m_arsenic_het" ,
                            ar_exposure_het= "abundance_ratio_log2_0_75_u_m_arsenic_nzo_0_75_u_m_arsenic_het" ,
                            p_ar_exposure = "abundance_ratio_p_value_0_75_u_m_arsenic_nzo_0_u_m_arsenic_nzo",
                            p_ar_0 = "abundance_ratio_p_value_0_u_m_arsenic_nzo_0_u_m_arsenic_het",
                            p_ar_exposure_het= "abundance_ratio_p_value_0_75_u_m_arsenic_nzo_0_75_u_m_arsenic_het") 
nod<- nod %>% dplyr::rename(ar_exposure= "abundance_ratio_log2_0_75_u_m_arsenic_nod_0_u_m_arsenic_nod",
                            ar_0 = "abundance_ratio_log2_0_u_m_arsenic_nod_0_u_m_arsenic_het" ,
                            ar_exposure_het= "abundance_ratio_log2_0_75_u_m_arsenic_nod_0_75_u_m_arsenic_het" ,
                            p_ar_exposure = "abundance_ratio_p_value_0_75_u_m_arsenic_nod_0_u_m_arsenic_nod",
                            p_ar_0 = "abundance_ratio_p_value_0_u_m_arsenic_nod_0_u_m_arsenic_het",
                            p_ar_exposure_het="abundance_ratio_p_value_0_75_u_m_arsenic_nod_0_75_u_m_arsenic_het" ) 
nod[,c(2:8)]<- lapply(nod[,c(2:8)], as.numeric) 
nod<-as.data.frame(nod)%>% dplyr::mutate(allele = "nod")
nzo[,c(2:8)]<- lapply(nzo[,c(2:8)], as.numeric) 
nzo<-as.data.frame(nzo)%>% dplyr::mutate(allele = "nzo")
both<-rbind(nod,nzo)


library(ggrepel)

gsea_theme2 <- theme(
  panel.grid.major = element_line(color = "gray70", size = 0.1), 
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  plot.title = element_text(hjust = 0.5), 
  #axis.text = element_text(size = 8),
  legend.key = element_rect(colour = "transparent", fill = "white"),
  axis.title = element_text(size = 16),
  #axis.title.y = element_blank(),
  axis.text.y=element_text(face = "bold",size = 14),
  axis.text.x=element_text(face = "bold",size = 14),
  strip.background = element_blank(),
  strip.text = element_text(color = "black",face = "bold",size = "14"))

#First set
nzo$diffexpressed <- "NO"
nzo$diffexpressed[nzo$ar_exposure > 0 & nzo$p_ar_exposure < 0.05] <- "UP"
nzo$diffexpressed[nzo$ar_exposure < 0 & nzo$p_ar_exposure < 0.05] <- "DOWN"


nod$diffexpressed <- "NO"
nod$diffexpressed[nod$ar_exposure > 0 & nod$p_ar_exposure < 0.05] <- "UP"
nod$diffexpressed[nod$ar_exposure < 0 & nod$p_ar_exposure < 0.05] <- "DOWN"

ggplot(nzo, aes(x= ar_exposure,y = -log10(p_ar_exposure), color = as.factor(diffexpressed),  label=gene_symbol))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NZO Exposure"))+ 
  scale_color_manual(values=c("#0072B2", "gray", "#D55E00")) +   geom_text_repel(min.segment.length = 0, seed = 42,segment.size = .2,max.overlaps = 20)
ggsave(
    filename = "~/Documents/arsenic_project/figures/nzo_exposure.pdf",
    width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)

ggplot(nod, aes(x= ar_exposure,y = -log10(p_ar_exposure), color = as.factor(diffexpressed),  label=gene_symbol))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NOD Exposure"))+ 
  scale_color_manual(values=c("#0072B2", "gray", "#D55E00")) +   geom_text_repel(min.segment.length = 0, seed = 42,segment.size = .2,max.overlaps = 20)
ggsave(
  filename = "~/Documents/arsenic_project/figures/nod_exposure.pdf",
  width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)


#Second set
nzo$diffexpressed <- "NO"
nzo$diffexpressed[nzo$ar_0 > 0 & nzo$p_ar_0 < 0.05] <- "UP"
nzo$diffexpressed[nzo$ar_0 < 0 & nzo$p_ar_0 < 0.05] <- "DOWN"


nod$diffexpressed <- "NO"
nod$diffexpressed[nod$ar_0 > 0 & nod$p_ar_0 < 0.05] <- "UP"
nod$diffexpressed[nod$ar_0 < 0 & nod$p_ar_0 < 0.05] <- "DOWN"

top_nzo<-  nzo %>% filter(ar_0 > 6) %>% pull(gene_symbol)
top_nzo<- top_nzo[top_nzo != "Txnrd1"]
top_nod<-  nod %>% filter(ar_0 > 6) %>% pull(gene_symbol)
top_nod<- top_nod[top_nod != "Txnrd1"]

ggplot(nzo, aes(x= ar_0,y = -log10(p_ar_0), color = as.factor(diffexpressed),  label=gene_symbol))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NZO Unexposed"))+ 
  scale_color_manual(values=c( "gray", "#D55E00")) +   geom_text_repel(min.segment.length = 0, seed = 42,segment.size = .2,max.overlaps = 20,hjust = 1,nudge_y = 1)
x_limits <- c( NA, 6.5)
ggsave(
  filename = "~/Documents/arsenic_project/figures/nzo_0_first.pdf",
  width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)

ggplot(nzo, aes(x= ar_0,y = -log10(p_ar_0), color = as.factor(diffexpressed),  label=ifelse(gene_symbol %in% top_nzo, gene_symbol, "")))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NZO Unexposed"))+ 
  scale_color_manual(values=c( "gray","#D55E00")) +    geom_text_repel(aes(
    segment.square = TRUE,
    segment.inflect = TRUE,
  ),
  force_pull   = 0, # do not pull toward data points
  nudge_y      = 0.05,
  nudge_x     = 0.05,
  direction    = "y",
  hjust        = 0,
  segment.size = 0.2,
  max.iter = 1e4, max.time = 1,max.overlaps = Inf,
  xlim = x_limits ,
  #segment.curvature = .1,
  box.padding = 1,
  segment.ncp = 3,
  segment.angle = 90,
  segment.curvature = .01,
  )
ggsave(
  filename = "~/Documents/arsenic_project/figures/nzo_0_second.pdf",
  width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)



ggplot(nod, aes(x= ar_0,y = -log10(p_ar_0), color = as.factor(diffexpressed),  label=gene_symbol))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NOD Unexposed"))+ 
  scale_color_manual(values=c( "gray", "#D55E00")) +   geom_text_repel(min.segment.length = 0, seed = 42,segment.size = .2,max.overlaps = 20)

ggsave(
  filename = "~/Documents/arsenic_project/figures/nod_0_first.pdf",
  width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)


ggplot(nod, aes(x= ar_0,y = -log10(p_ar_0), color = as.factor(diffexpressed),  label=ifelse(gene_symbol %in% top_nod, gene_symbol, "")))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NOD Unexposed"))+ 
  scale_color_manual(values=c( "gray","#D55E00")) +    geom_text_repel(aes(
    segment.square = TRUE,
    segment.inflect = TRUE,
  ),
  force_pull   = 0, # do not pull toward data points
  nudge_y      = 0.05,
  nudge_x     = 0.05,
  direction    = "y",
  hjust        = 0,
  segment.size = 0.2,
  max.iter = 1e4, max.time = 1,max.overlaps = Inf,
  xlim = x_limits ,
  #segment.curvature = .1,
  box.padding = 1,
  segment.ncp = 3,
  segment.angle = 90,
  segment.curvature = .01,
  )
ggsave(
  filename = "~/Documents/arsenic_project/figures/nod_0_second.pdf",
  width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)





#Third set
nzo$diffexpressed <- "NO"
nzo$diffexpressed[nzo$ar_exposure_het > 0 & nzo$p_ar_exposure_het < 0.05] <- "UP"
nzo$diffexpressed[nzo$ar_exposure_het < 0 & nzo$p_ar_exposure_het < 0.05] <- "DOWN"


nod$diffexpressed <- "NO"
nod$diffexpressed[nod$ar_exposure_het > 0 & nod$p_ar_exposure_het < 0.05] <- "UP"
nod$diffexpressed[nod$ar_exposure_het < 0 & nod$p_ar_exposure_het < 0.05] <- "DOWN"

top_nzo<-  nzo %>% filter(ar_exposure_het > 6) %>% pull(gene_symbol)
top_nzo<- top_nzo[top_nzo != "Txnrd1"]
top_nod<-  nod %>% filter(ar_exposure_het > 6) %>% pull(gene_symbol)
top_nod<- top_nod[top_nod != "Txnrd1"]


ggplot(nzo, aes(x= ar_exposure_het,y = -log10(p_ar_exposure_het), color = as.factor(diffexpressed),  label=gene_symbol))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NZO/Het Exposure"))+ 
  scale_color_manual(values=c( "gray", "#D55E00")) +   geom_text_repel(min.segment.length = 0, seed = 42,segment.size = .2,max.overlaps = 20)
ggsave(
  filename = "~/Documents/arsenic_project/figures/nzo_het_first.pdf",
  width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)

ggplot(nzo, aes(x= ar_exposure_het,y = -log10(p_ar_exposure_het), color = as.factor(diffexpressed),  label=ifelse(gene_symbol %in% top_nzo, gene_symbol, "")))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NZO/Het Exposure"))+ 
  scale_color_manual(values=c("gray", "#D55E00")) +    geom_text_repel(aes(
    segment.square = TRUE,
    segment.inflect = TRUE,
  ),
  force_pull   = 0, # do not pull toward data points
  nudge_y      = 0.05,
  nudge_x     = 0.05,
  direction    = "y",
  hjust        = 0,
  segment.size = 0.2,
  max.iter = 1e4, max.time = 1,max.overlaps = Inf,
  xlim = x_limits ,
  #segment.curvature = .1,
  box.padding = 1,
  segment.ncp = 3,
  segment.angle = 90,
  segment.curvature = .01,
  )
ggsave(
  filename = "~/Documents/arsenic_project/figures/nzo_het_second.pdf",
  width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)



ggplot(nod, aes(x= ar_exposure_het,y = -log10(p_ar_exposure_het), color = as.factor(diffexpressed),  label=gene_symbol))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NOD/Het Exposure"))+ 
  scale_color_manual(values=c( "gray", "#D55E00")) +   geom_text_repel(min.segment.length = 0, seed = 42,segment.size = .2,max.overlaps = 20)
ggsave(
  filename = "~/Documents/arsenic_project/figures/nod_het_first.pdf",
  width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)

ggplot(nod, aes(x= ar_exposure_het,y = -log10(p_ar_exposure_het), color = as.factor(diffexpressed),  label=ifelse(gene_symbol %in% top_nod, gene_symbol, "")))+geom_point(alpha = .9,shape = 19, stroke=NA)+theme(legend.position = "none")+gsea_theme2+
  gsea_theme2+labs(color = "P value")+ geom_vline(xintercept=c(0,0), col="darkgray",lty = 2,lwd = 1) +geom_hline(yintercept=-log10(0.05), col="darkgray",lty = 2,lwd = 1)+labs(y=expression(-log[10]~"P"),x = expression(log[2]~"Abundance Ratio NOD/Het Exposure"))+ 
  scale_color_manual(values=c( "gray","#D55E00")) +    geom_text_repel(aes(
    segment.square = TRUE,
    segment.inflect = TRUE,
  ),
  force_pull   = 0, # do not pull toward data points
  nudge_y      = 0.05,
  nudge_x     = 0.05,
  direction    = "y",
  hjust        = 0,
  segment.size = 0.2,
  max.iter = 1e4, max.time = 1,max.overlaps = Inf,
  xlim = x_limits ,
  #segment.curvature = .1,
  box.padding = 1,
  segment.ncp = 3,
  segment.angle = 90,
  segment.curvature = .01,
  )
ggsave(
  filename = "~/Documents/arsenic_project/figures/nod_het_second.pdf",
  width = 7, height = 7, dpi = 300, units = "in", device='pdf',useDingbats = F
)


