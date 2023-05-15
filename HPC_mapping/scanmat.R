library(tidyverse)
library(purrr)
tbl <-
    list.files(path = "/pod/2/reinholdt-lab/Callan/spring_2023/output",
               pattern = "*.rds",
               full.names = T)
tbl2<-
    map_dfc(tbl, ~readRDS(.), full.names = T)
cnames<- (as.data.frame(tbl) %>% mutate(feature = unglue::unglue_vec(tbl,"{}output/{x}.{}")))
colnames(tbl2)<- cnames$feature

saveRDS(tbl2,"/pod/2/reinholdt-lab/Callan/spring_2023/output2/scanmatrix.rds")
