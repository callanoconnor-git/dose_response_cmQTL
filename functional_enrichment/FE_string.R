# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("STRINGdb")

library(STRINGdb)
stringdat<- read_table("~/Downloads/10090.protein.links.v11.5 (1).txt")
dim(stringdat)
#https://www.youtube.com/watch?v=JEaNhUq_8rA


#define response cmQTL
getdir<- readRDS("~/Documents/arsenic_project/getdirection.rds")
getdir_base <- gs_sum3 %>% filter(!trait %in% getdir$trait) %>% dplyr::mutate(par = unglue_vec(trait,"{x}_{}"))
getdir_base<- getdir_base %>% filter(par == "c" | par == "d")
base_exp_qtl<- res_plot %>% filter(gene %in% getdir_base$ensembl_id)

gs_sum3<- readRDS("~/Documents/arsenic_project/data/geneset_ctd_candidates.rds")

responsefromqtl<- gs_sum3 %>% filter(!trait %in% getdir_base$trait)
response_set<- responsefromqtl %>% dplyr::select(ensembl_id,Frequency ) %>% unique()
response_set<- res_plot %>% dplyr::filter(gene %in% unique(responsefromqtl$ensembl_id))
#How many response cmQTL
length(responsefromqtl$trait)
#How many unique genes
length(unique(response_set$gene))

###
#Run PPI/STRING
string_db <- STRINGdb$new( version="11.5", species=10090,
                           score_threshold=700, input_directory="")
response_mapped <- string_db$map( response_set, "gene", removeUnmappedRows = TRUE )
hits_response <- response_mapped$STRING_id
#write.csv(response_mapped,"~/Documents/arsenic_project/results/FE_cmQTL.csv")
#plota<- string_db$plot_network( hits_response )
enrichment<- string_db$get_enrichment(response_mapped$STRING_id)
#write.csv(enrichment,"~/Documents/arsenic_project/results/cmQTL_enrichment.csv")

clusterList<- string_db$get_clusters(response_mapped$STRING_id, algorithm = 'fastgreedy')
min_len<- sapply(clusterList,function(x)length(x)>1)
clusterList<- clusterList[min_len]
#string_db$plot_network(clusterList[[1]])

ggplot(subset(enrichment, category == 'MPO'),aes(number_of_genes,reorder(description,description),fill = p_value))+
  geom_bar(stat = "identity")+
geom_text(aes(label = number_of_genes),color = "black",hjust = -0.5)+
  theme_bw()


#compartments<- subset(enrichment, category == "COMPARTMENTS")%>% dplyr::mutate(genese="GOCC") %>% slice_min(fdr,n=12)
process<- subset(enrichment, category == "Process")%>% dplyr::mutate(genese="GOBP")%>% slice_min(fdr,n=12)
component<- subset(enrichment, category == "Component")%>% dplyr::mutate(genese="GOCC")%>% slice_min(fdr,n=12)
kegg<- subset(enrichment, category == "KEGG")%>% dplyr::mutate(genese="KEGG")%>% slice_min(fdr,n=12)
wiki<- subset(enrichment, category == "WikiPathways")%>% dplyr::mutate(genese="WP")%>% slice_min(fdr,n=12)
rctm<- subset(enrichment, category == "RCTM") %>% dplyr::mutate(genese="REACTOME")%>% slice_min(fdr,n=12)
go_mf<- subset(enrichment, category == "Function") %>% dplyr::mutate(genese="GOMF")%>% slice_min(fdr,n=12)
tissue<- subset(enrichment, category == "TISSUES") %>% dplyr::mutate(genese="Tissue")%>% slice_min(fdr,n=12)
#pmid<- subset(enrichment, category == "PMID") %>% dplyr::mutate(genese="PMID")%>% slice_min(fdr,n=12)
mpo<- subset(enrichment, category == "MPO") %>% dplyr::mutate(genese="MPO")%>% slice_min(fdr,n=12)
allbind<- rbind(process,component,kegg,wiki,tissue,rctm,go_mf,mpo)%>% slice_min(fdr,n=12)
#write.csv(allbind,"~/Documents/arsenic_project/figures/allpeaks_string/allbind.csv")
#write.csv(enrichment,"~/Documents/arsenic_project/figures/allpeaks_string/cmQTL_candidate_enrichment.csv")

ggplot(subset(enrichment, category == 'Process'),aes(number_of_genes,description,fill = p_value))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = number_of_genes),color = "black",hjust = -0.5)+
  theme_bw()

library(ggtext)
library(glue)
allbind<- allbind %>% dplyr::group_by(category) %>% arrange(category,fdr) %>%  ungroup()%>%
         dplyr::mutate(
           color = c(rep("#000000",12),
                    rep("#E69F00",12),
                    rep("#56B4E9",12),
                    rep("#009E73",13),
                    rep("#F0E442",12),
                    rep("#0072B2",12),
                    rep("#D55E00",12),
                    rep("#CC79A7",12)))
allbind$color<- as.factor(allbind$color)
allbind$description<- as.factor(allbind$description)
allbind$description2<- paste0("<i style=\'color:", allbind$color, "\'>", allbind$description, "</i>")
plot_1<- allbind %>%
  ggplot(   aes(number_of_genes,description2,fill = fdr))+
  geom_bar(stat = "identity")+
  #geom_text(aes(label = number_of_genes),color = "black",hjust = -0.5)+theme_bw()+
  theme(axis.text.x  = element_markdown(size = 6,angle = 90,hjust=0.95,vjust=0.2))+xlab(~"")+xlab(~"Number of Genes")+coord_flip()
#ggsave(filename = "~/Documents/arsenic_project/figures/allpeaks_string/functional_enrichment.pdf",plot_1, width = 8, height = 4, dpi = 300, units = "in", device='pdf',useDingbats = F)

###
subnet<- string_db$get_subnetwork(hits_response)
all_edges<- string_db$get_interactions(response_mapped$STRING_id)
neighbors<- string_db$get_neighbors(response_mapped$STRING_id)
neighbors2<-left_join(neighbors %>% as.data.frame(),response_mapped, by = c("." = "STRING_id"))
go1<- left_join(all_edges,response_mapped, by = c("from" = "STRING_id"))
go1<- left_join(go1,response_mapped[,-c(2:8)], by = c("to" = "STRING_id"))
go1<- go1 %>% dplyr::rename(fromgene = "gene.x",
                            togene = "gene.y")
allgenes<- read_delim("~/Downloads/MGI_Gene_Model_Coord.rpt.txt")
allgenes<- allgenes[,c(3,11)]
colnames(allgenes)<- c("symbol","gene")
go2<- left_join(go1,allgenes, by = c("fromgene" = "gene" )) 
go2<- left_join(go2,allgenes, by = c("togene" = "gene" )) 
go2<- go2 %>% dplyr::select(symbol.x,symbol.y,combined_score) %>% dplyr::rename(from = "symbol.x",
                            to = "symbol.y")
golist<- c(go2$from,go2$to) %>% unique()
go2<- go2 %>% unique()
library(igraph)


g<- graph_from_data_frame(go2,directed=FALSE,vertices = golist)
plot(g)
# see how many proteins do you have    
vcount(g)
degree(g)
sort(degree(g),decreasing = T)

#https://stackoverflow.com/questions/26218900/making-igraph-clearer-to-read
pdf("~/Documents/arsenic_project/figures/igraph1_test.pdf", width = 8, height = 11)
l <- layout.fruchterman.reingold(g, niter=5000, area=vcount(g)^4*10)

plot(g, layout=l, 
     edge.arrow.size=0.5, 
     vertex.label.cex=0.75, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     vertex.la=2,
     vertex.shape="circle", 
     vertex.size=10, 
     vertex.label.color="black", 
     edge.width=1)
dev.off()

# find top 200 proteins with the highest degree
top.degree.verticies <- names(tail(sort(degree(full.graph)), 200))
# count the number of proteins in it
vcount(full.graph)



# gsea_res<- readRDS("~/Documents/arsenic_project/data/gsea_dose_scores.rds")
# gsea_res<- gsea_res%>% dplyr::mutate(geneset = unglue::unglue_vec(pathway,"{x}_{}")) %>% dplyr::mutate(path = gsub( "_"," ", unglue::unglue_vec(pathway,"{}_{x}") ))
# 

go_path<- subset(enrichment, category == "MPO")
# go_path<- subset(enrichment,  category == "Process" & description =="Binding" )
go_path<- subset(enrichment,  category == "MPO" & description =="Abnormal cell death" )

string_id<- str_split(go_path$inputGenes,pattern = ',')
genesymbol<- str_split(go_path$preferredNames,pattern = ',')
go_path<- data.frame(nodes = string_id,genesymbol = genesymbol)
colnames(go_path)<-c('nodes','genesymbol')
#string_db$plot_network(go_path$nodes)

all_edges<- string_db$get_interactions(go_path$nodes)
all_edges<-left_join(all_edges,go_path,by=c('from'='nodes'))
all_edges<-left_join(all_edges,go_path,by=c('to'='nodes'))
all_edges<-unique(all_edges)

all_edges<-all_edges[,4:5]
colnames(all_edges)<-c('from','to')

library(igraph)
g<- graph_from_data_frame(all_edges,directed=FALSE,vertices = go_path$genesymbol)
g
vcount(g)
degree(g)
sort(degree(g),decreasing = T)
E(g)
V(g)$color
plot(g)

#######
nod_ip<- (nod %>% filter(p_ar_exposure < .05 & ar_exposure > .75))
string_db <- STRINGdb$new( version="11.5", species=10090,
                           score_threshold=400, input_directory="")
nod_mapped <- string_db$map( nod_ip, "gene_symbol", removeUnmappedRows = TRUE )
nod_response <- nod_mapped$STRING_id
#write.csv(nod_mapped,"~/Documents/arsenic_project/results/FE_NOD_IP.csv")
#plota<- string_db$plot_network( nod_response )

enrichment<- string_db$get_enrichment(nod_mapped$STRING_id)

ggplot(subset(enrichment, category == 'Component'),aes(number_of_genes,description,fill = p_value))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = number_of_genes),color = "black",hjust = -0.5)+
  theme_bw()

all_edges<- string_db$get_interactions(nod_mapped$STRING_id)
neighbors<- string_db$get_neighbors(nod_mapped$STRING_id)
neighbors2<-left_join(neighbors %>% as.data.frame(),nod_mapped, by = c("." = "STRING_id"))
go1<- left_join(all_edges,nod_mapped, by = c("from" = "STRING_id"))
go1<- left_join(go1,nod_mapped[,-c(2:18)], by = c("to" = "STRING_id"))
go1<- go1 %>% dplyr::rename(fromgene = "gene_symbol.x",
                            togene = "gene_symbol.y")
allgenes<- read_delim("~/Downloads/MGI_Gene_Model_Coord.rpt.txt")
allgenes<- allgenes[,c(3,11)]
colnames(allgenes)<- c("symbol","gene")
go2<- left_join(go1,allgenes, by = c("fromgene" = "gene" )) 
go2<- left_join(go2,allgenes, by = c("togene" = "gene" )) 
go2<- go2 %>% dplyr::select(fromgene,togene,combined_score) %>% dplyr::rename(from = "fromgene",
                                                                                to = "togene")
golist<- c(go2$from,go2$to) %>% unique()
go2<- go2 %>% unique()
library(igraph)


g<- graph_from_data_frame(go2,directed=FALSE,vertices = golist)
#plot(g)
# see how many proteins do you have    
vcount(g)
degree(g)
sort(degree(g),decreasing = T)

#https://stackoverflow.com/questions/26218900/making-igraph-clearer-to-read
pdf("~/Documents/arsenic_project/figures/igraph1_NOD.pdf", width = 8, height = 11)
l <- layout.fruchterman.reingold(g, niter=5000, area=vcount(g)^4*10)

plot(g, layout=l, 
     edge.arrow.size=0.5, 
     vertex.label.cex=0.75, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     vertex.la=2,
     vertex.shape="circle", 
     vertex.size=10, 
     vertex.label.color="black", 
     edge.width=1)
dev.off()
#######
#######
nzo_ip<- (nzo %>% filter(p_ar_exposure < .05 & ar_exposure > .75))
string_db <- STRINGdb$new( version="11.5", species=10090,
                           score_threshold=400, input_directory="")
nzo_mapped <- string_db$map( nzo_ip, "gene_symbol", removeUnmappedRows = TRUE )
nzo_response <- nzo_mapped$STRING_id
write.csv(nzo_mapped,"~/Documents/arsenic_project/results/FE_NZO_IP.csv")
plota<- string_db$plot_network( nzo_response )

enrichment<- string_db$get_enrichment((nzo_mapped$STRING_id))
#write.csv(enrichment,"~/Documents/arsenic_project/results/FE_NZO_IP_enrichment.csv")

ggplot(subset(enrichment, category == 'Component'),aes(number_of_genes,description,fill = p_value))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = number_of_genes),color = "black",hjust = -0.5)+
  theme_bw()

go_path<- subset(enrichment, description = "Process")
go_path<- subset(enrichment, description =="Protein binding" & category == "Function")

string_id<- str_split(go_path$inputGenes,pattern = ',')
genesymbol<- str_split(go_path$preferredNames,pattern = ',')
go_path<- data.frame(nodes = string_id,genesymbol = genesymbol)
colnames(go_path)<-c('nodes','genesymbol')
string_db$plot_network(go_path$nodes)

all_edges<- string_db$get_interactions(go_path$nodes)
all_edges<-left_join(all_edges,go_path,by=c('from'='nodes'))
all_edges<-left_join(all_edges,go_path,by=c('to'='nodes'))
all_edges<-unique(all_edges)

all_edges<-all_edges[,4:5]
colnames(all_edges)<-c('from','to')

library(igraph)
g<- graph_from_data_frame(all_edges,directed=FALSE,vertices = go_path$genesymbol)
g
vcount(g)
degree(g)
sort(degree(g),decreasing = T)
E(g)
V(g)$color
plot(g)


#######

allele_2<- readRDS("~/Documents/arsenic_project/data/allele_2.rds")
allele_de<- as.data.frame(allele_2) %>% filter(padj < .05) %>% dplyr::mutate(gene = rownames(.))
string_db <- STRINGdb$new( version="11.5", species=10090,
                           score_threshold=400, input_directory="")
allele_mapped <- string_db$map( allele_de, "gene", removeUnmappedRows = TRUE )
allele_response <- allele_mapped$STRING_id
#plota<- string_db$plot_network( hits_response )
#write.csv(allele_mapped,"~/Documents/arsenic_project/results/FE_allele_DE.csv")
enrichment<- string_db$get_enrichment((response_mapped$STRING_id))

ggplot(subset(enrichment, category == 'RCTM'),aes(number_of_genes,description,fill = p_value))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = number_of_genes),color = "black",hjust = -0.5)+
  theme_bw()

go_path<- subset(enrichment, description = "Process")
go_path<- subset(enrichment, description =="Protein binding" & category == "Function")

all_edges<- string_db$get_interactions(allele_mapped$STRING_id)
neighbors<- string_db$get_neighbors(allele_mapped$STRING_id)
neighbors2<-left_join(neighbors %>% as.data.frame(),allele_mapped, by = c("." = "STRING_id"))
go1<- left_join(all_edges,allele_mapped, by = c("from" = "STRING_id"))
go1<- left_join(go1,allele_mapped[,], by = c("to" = "STRING_id"))
go1<- go1 %>% dplyr::rename(fromgene = "gene.x",
                            togene = "gene.y")
allgenes<- read_delim("~/Downloads/MGI_Gene_Model_Coord.rpt.txt")
allgenes<- allgenes[,c(3,11)]
colnames(allgenes)<- c("symbol","gene")
go2<- left_join(go1,allgenes, by = c("fromgene" = "gene" )) 
go2<- left_join(go2,allgenes, by = c("togene" = "gene" )) 
go2<- go2 %>% dplyr::select(symbol.x,symbol.y,combined_score) %>% dplyr::rename(from = "symbol.x",
                                                                                to = "symbol.y")
golist<- c(go2$from,go2$to) %>% unique()
go2<- go2 %>% unique()
library(igraph)


g<- graph_from_data_frame(go2,directed=FALSE,vertices = golist)
#plot(g)
# see how many proteins do you have    
vcount(g)
degree(g)
sort(degree(g),decreasing = T)

#https://stackoverflow.com/questions/26218900/making-igraph-clearer-to-read
pdf("~/Documents/arsenic_project/figures/igraph1_allele2.pdf", width = 8, height = 11)
l <- layout.fruchterman.reingold(g, niter=5000, area=vcount(g)^4*10)

plot(g, layout=l, 
     edge.arrow.size=0.5, 
     vertex.label.cex=0.75, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     vertex.la=2,
     vertex.shape="circle", 
     vertex.size=10, 
     vertex.label.color="black", 
     edge.width=1)
dev.off()


# 
# 
# baseline<- gs_sum2 %>% filter(trait %in% getdir_base$trait)
# baseline_set<- res_plot %>% dplyr::filter(gene %in% unique(baseline$ensembl_id))
# string_db <- STRINGdb$new( version="11.5", species=10090,
#                            score_threshold=700, input_directory="")
# base_mapped <- string_db$map( baseline_set, "gene", removeUnmappedRows = TRUE )
# base_response<- base_mapped$STRING_id
# 
# enrichment<- string_db$get_enrichment(base_mapped$STRING_id)
# 
# #write.csv(enrichment,"process<- subset(enrichment, category == "Process")%>% dplyr::mutate(genese="GOBP")%>% slice_min(fdr,n=12)
# #compartments<- subset(enrichment, category == "COMPARTMENTS")%>% dplyr::mutate(genese="GOCC") %>% slice_min(fdr,n=12)
# component<- subset(enrichment, category == "Component")%>% dplyr::mutate(genese="GOCC")%>% slice_min(fdr,n=12)
# process<- subset(enrichment, category == "Process")%>% dplyr::mutate(genese="GOBP")%>% slice_min(fdr,n=12)
# kegg<- subset(enrichment, category == "KEGG")%>% dplyr::mutate(genese="KEGG")%>% slice_min(fdr,n=12)
# wiki<- subset(enrichment, category == "WikiPathways")%>% dplyr::mutate(genese="WP")%>% slice_min(fdr,n=12)
# rctm<- subset(enrichment, category == "RCTM") %>% dplyr::mutate(genese="REACTOME")%>% slice_min(fdr,n=12)
# go_mf<- subset(enrichment, category == "Function") %>% dplyr::mutate(genese="GOMF")%>% slice_min(fdr,n=12)
# tissue<- subset(enrichment, category == "TISSUES") %>% dplyr::mutate(genese="Tissue")%>% slice_min(fdr,n=12)
# #pmid<- subset(enrichment, category == "PMID") %>% dplyr::mutate(genese="PMID")%>% slice_min(fdr,n=12)
# mpo<- subset(enrichment, category == "MPO") %>% dplyr::mutate(genese="MPO") 
# allbind<- rbind(process,component,kegg,wiki,tissue,rctm,go_mf,mpo)
# #write.csv(allbind,"~/Documents/arsenic_project/figures/allpeaks_string/allbind.csv")
# #write.csv(enrichment,"~/Documents/arsenic_project/figures/allpeaks_string/cmQTL_candidate_enrichment.csv")
# component<- subset(enrichment, category == "Component")%>% dplyr::mutate(genese="GOCC") 
# process<- subset(enrichment, category == "Process")%>% dplyr::mutate(genese="GOBP") 
# kegg<- subset(enrichment, category == "KEGG")%>% dplyr::mutate(genese="KEGG") 
# wiki<- subset(enrichment, category == "WikiPathways")%>% dplyr::mutate(genese="WP") 
# rctm<- subset(enrichment, category == "RCTM") %>% dplyr::mutate(genese="REACTOME") 
# go_mf<- subset(enrichment, category == "Function") %>% dplyr::mutate(genese="GOMF") 
# tissue<- subset(enrichment, category == "TISSUES") %>% dplyr::mutate(genese="Tissue") 
# #pmid<- subset(enrichment, category == "PMID") %>% dplyr::mutate(genese="PMID") 
# mpo<- subset(enrichment, category == "MPO") %>% dplyr::mutate(genese="MPO") 
# allbind<- rbind(process,component,kegg,wiki,tissue,rctm,go_mf,mpo)
# 
# ggplot(subset(enrichment, category == 'Process'),aes(number_of_genes,description,fill = p_value))+
#   geom_bar(stat = "identity")+
#   geom_text(aes(label = number_of_genes),color = "black",hjust = -0.5)+
#   theme_bw()
# 
# library(ggtext)
# library(glue)
# allbind<- allbind %>% dplyr::group_by(category) %>% arrange(category,fdr) %>%  ungroup()%>%
#   dplyr::mutate(
#     color = c(rep("#000000",13),
#               rep("#E69F00",18),
#              # rep("#56B4E9",12),
#              # rep("#009E73",12),
#               rep("#F0E442",5),
#               #rep("#0072B2",12),
#               rep("#D55E00",1),
#               rep("#CC79A7",8)))
# allbind$color<- as.factor(allbind$color)
# allbind$description<- as.factor(allbind$description)
# allbind$description2<- paste0("<i style=\'color:", allbind$color, "\'>", allbind$description, "</i>")
# plot_1<- allbind %>%
#   ggplot(   aes(number_of_genes,description2,fill = fdr))+
#   geom_bar(stat = "identity")+
#   #geom_text(aes(label = number_of_genes),color = "black",hjust = -0.5)+theme_bw()+
#   theme(axis.text.x  = element_markdown(size = 6,angle = 90,hjust=0.95,vjust=0.2))+xlab(~"")+xlab(~"Number of Genes")+coord_flip()
# 
# 
# 
# 
# 
