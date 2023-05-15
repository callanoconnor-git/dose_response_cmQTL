library(tidyverse)
library(lme4)


dat<- readr::read_csv("~/Desktop/relativeabundance/genotypeeffect.csv")%>% janitor::clean_names()
dat<- dat %>% filter(sex != "other" & animal_genotype != "WT"& animal_genotype != "MT") 
dat<- dat %>%  dplyr::mutate(conc = unglue::unglue_vec(sample_name, "{}_{}_{x}_{}"))
dat<- dat %>%  dplyr::mutate(rep = unglue::unglue_vec(sample_name, "{}_{}_{}_{x}"))
dat<- dat %>%  dplyr::mutate(ind = unglue::unglue_vec(sample_name, "{}_{x}_{}_{}"))

dat<- dat %>% filter(ind != 4543)
dat<- dat %>% dplyr::mutate(indcolor = ifelse(animal_genotype == "NOD", "#1111FF",
                                              ifelse(animal_genotype == "NZO", "#00CCFF",
                                                     ifelse(animal_genotype == "HET", "gray","gray"
                                                     ))))
indcolor<- dat$indcolor

txnrd1_theme <- theme(
  legend.position = "none",
  #  axis.text.y=element_text(size=6,face = "bold"),
  text = element_text(face = "bold"),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  plot.title = element_text(hjust = 0.5), 
  axis.text = element_text(size = 20),
  axis.title = element_text(size = 20),
  legend.key.height= unit(1, 'cm'),
  legend.key.width= unit(1.5, 'cm'),
  legend.text = element_text(size=20),
  legend.title = element_text(size = 20))

mod1<- lm(area_ratio_to_mean ~ animal_genotype+sex, data = dat %>% filter(animal_genotype != "HET" & conc == 0 & ind != 4543)  )
anova(mod1)

mod2<- lm(area_ratio_to_mean ~ conc+animal_genotype+sex + conc:animal_genotype+conc:sex, data = dat %>% filter(animal_genotype != "HET" & ind != 4543  ))
anova(mod2)

mod3<- lm(area_ratio_to_mean ~ animal_genotype+sex, data = dat %>% filter(animal_genotype != "HET" & conc !=0 & ind != 4543)  )
anova(mod3)

dat$conc[dat$conc =='0-75']<- 0.75
dat$conc<- as.numeric(dat$conc)


datforperms<- dat %>% filter(animal_genotype != "HET" & conc == 0 & ind != 4543) %>% dplyr::select(animal_genotype,area_ratio_to_mean)
set.seed(101) ## for reproducibility
nsim <- 1000
res <- numeric(nsim) ## set aside space for results
for (i in 1:nsim) {
  ## standard approach: scramble response value
  perm <- sample(nrow(datforperms))
  bdat <- transform(datforperms,area_ratio_to_mean=area_ratio_to_mean[perm])
  ## compute & store difference in means; store the value
  mod1<- lm(area_ratio_to_mean ~ animal_genotype, data = bdat )
  an1<- anova(mod1)
  res[i] <-  an1$`Pr(>F)`[1]
}

##Callans model
mod1<- lm(area_ratio_to_mean ~ animal_genotype, data = dat %>% filter(animal_genotype != "HET" & conc == 0 & ind != 4543)  )
callan<- anova(mod1)

mean(res) > callan$`Pr(>F)`[1]


mod1<- lm(area_ratio_to_mean ~ animal_genotype, data = dat %>% filter(animal_genotype != "HET" & conc == 0 & ind != 4543)  )
anova(mod1)


obs <- mean(datforperms$area_ratio_to_mean[datforperms$animal_genotype=="NOD"])-
  mean(datforperms$area_ratio_to_mean[datforperms$animal_genotype=="NZO"])
## append the observed value to the list of results
res <- c(res,obs)

res

hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")

mean(abs(res)>=abs(obs)) 
############


dat %>% filter( conc == 0) %>%
ggplot(aes(x = as.factor(animal_genotype), y = area_ratio_to_mean, fill = as.factor(animal_genotype)))+
  geom_boxplot(alpha = .7,outlier.shape = NA)+geom_jitter(aes(col = as.factor(animal_genotype),shape = as.factor(sex)),size = 3,width = .12)+
  scale_color_manual(values = c("gray","#0072B2","#56B4E9"))+
  scale_fill_manual(values = c("gray","#0072B2","#56B4E9"))+txnrd1_theme+
  xlab(~"Haplotype")+labs(color = "Haplotype", fill = "Haplotype", shape = "Sex")+
  ylab(~"TXNRD1 Levels 0 uM")

dat %>% filter( conc == 0.75) %>%
  ggplot(aes(x = as.factor(animal_genotype), y = area_ratio_to_mean, fill = as.factor(animal_genotype)))+
  geom_boxplot(alpha = .7,outlier.shape = NA)+geom_jitter(aes(col = as.factor(animal_genotype),shape = as.factor(sex)),size = 3,width = .12)+
  scale_color_manual(values = c("gray","#0072B2","#56B4E9"))+
  scale_fill_manual(values = c("gray","#0072B2","#56B4E9"))+txnrd1_theme+
  xlab(~"Haplotype")+labs(color = "Haplotype", fill = "Haplotype", shape = "Sex")+
  ylab(~"TXNRD1 Levels 0.75 uM")

dat %>%
  ggplot(aes(x = interaction(animal_genotype,conc,sep = "!"), y = area_ratio_to_mean, fill = as.factor(animal_genotype)))+
  geom_boxplot(alpha = .7,outlier.shape = NA)+geom_jitter(aes(x = interaction(animal_genotype,conc, sep = "!"),col = as.factor(animal_genotype),shape = as.factor(sex)),size = 3,width = .12)+
  scale_color_manual(values = c("gray","#0072B2","#56B4E9"))+
  scale_fill_manual(values = c("gray","#0072B2","#56B4E9"))+
  txnrd1_theme+
  xlab(~"Haplotype")+labs(color = "Haplotype", fill = "Haplotype", shape = "Sex")+
  ylab(~"Relative Abundance") + scale_x_discrete(guide = ggh4x::guide_axis_nested(delim = "!"), name = "Concentration")
#ggsave(filename = "~/Documents/arsenic_project/figures/txnrd1_rel_adbundance.pdf",width = 5, height = 4, dpi = 300, units = "in", device='pdf',useDingbats = F)

