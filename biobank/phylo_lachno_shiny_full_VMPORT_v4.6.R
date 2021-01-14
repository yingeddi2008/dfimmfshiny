rm(list=ls())
library(ape)
library(tidyverse)
library(RPostgreSQL)
library(ggtree)
library(shiny)
library(stringr)
library(yingtools2)
library(shinythemes)
library(ggtree)
library(grid)
library(gridExtra)
library(DT)


# original get tables --------------------------------------------------------------------

#load("/Volumes/pamer-lab/Eric.Littmann/Projects/msk_isolate_library/make_tree/lachno_tree_int.RData")
if(F){
killDbConnections <- function (){
  
  all_cons <- dbListConnections(dbDriver("PostgreSQL"),host="iski0035",dbname="isolate_library",
                                user="postgres",password="pamerlab1625")
  
  print(all_cons)
  
  for(con in all_cons)
    +  dbDisconnect(con)
  
  print(paste(length(all_cons), " connections killed."))
  
}
#killDbConnections()

#updated tree
tree <- read.tree("/Volumes/pamer-lab/Eric.Littmann/Projects/msk_isolate_library/shiny_things/refs/full_open_biome_msk_phyml_tree.txt")

con <- dbConnect(dbDriver("PostgreSQL"),
                  host="128.135.41.183",
                # host="10.151.15.23",
                 dbname="dfi_commensal_library",
                 user="ericlittmann",
                 password="dfibugs")
#dbListTables(con)

#get that isolate info
exclude <- read.csv(file="/Volumes/pamer-lab/Eric.Littmann/Projects/lachno_rebuttal/seq_id_exclude_duplicates.csv",
                    stringsAsFactors = F)
lookup <- tbl(con,"lookup_table") %>% collect()
# matt_tax <- tbl(con,"matt_tax") %>% collect() %>%
#    left_join(lookup) %>%
#    filter(is_lachno==1,
#           seq_id %!in% exclude$seq_id_exlcude)
blast_raw <-tbl(con,"wgs_16s_patric_blast_results_copynum")

blast_top <- blast_raw %>%
  filter(length>0) %>%
  collect() %>%
  group_by(seq_id,which_copy_num) %>%
  arrange(desc(bitscore)) %>%
  dplyr::slice(1:5) %>%
  ungroup() %>%
  group_by(seq_id) %>%
  dplyr::select(seq_id,query_acc,
         which_copy_num,rank,phylum,
         class,order,family,genus,
         species,bitscore,pident,length) %>%
  group_by(seq_id) %>%
  mutate(total_copy=max(which_copy_num)) %>%
  group_by(query_acc) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  group_by(seq_id) %>%
  filter(all(phylum==unique(phylum)),
         all(genus==unique(genus))) %>%
  filter(which_copy_num==1) %>%
  dplyr::slice(1) %>%
  ungroup()


blast <- tbl(con,"clean_taxonomy_lookup") %>% 
  collect() %>%
  left_join(blast_top %>%
              dplyr::select(seq_id,length))
  
#dbListTables(con)
 seq <- tbl(con,"prokka_sequences") %>% collect()
# 
# seqs <- tbl(con,"prokka_annotations") %>%
#   filter(seq_id %in% !!blast$seq_id) %>%
#   left_join(seq) %>%
#   collect()
# 

meta <- tree$tip.label %>%
  as.data.frame() %>%
  mutate(seq_id=gsub("\\_[mo].+.+","",.),
         length=str_extract(pattern="[0-9]+$",.),
         cohort=str_extract(.,pattern="msk|open_biome")) %>%
  left_join(lookup) %>%
  select(-length) %>%
  left_join(blast) %>%
  mutate(family=ifelse(grepl(species,pattern="Ruminoc"),"Ruminococcaceae",family),
         int_id=gsub("MSK\\.","Msk",msk_id),
         int_id=gsub("\\.","_",int_id),
         wiki_link=paste0("http://128.135.41.103/wiki/index.php/",int_id)) %>%
  select(-int_id) %>%
  mutate(both_ids=paste(msk_id,species))

colnames(meta)[1] <- "tiplab"


#make big full taxonomy table
cleantax <- tbl(con,"clean_taxonomy_lookup_v2") %>% collect() %>%
  select(-row.names)

blast_int <- blast_top %>%
  filter(seq_id %!in% cleantax$seq_id) %>%
  select(seq_id,phylum,class,order,family,genus,species) %>%
  left_join(lookup %>% select(seq_id,msk_id,donor_id)) %>%
  mutate(mixed=1,ambiguous=1,comment="genome may be contaminated")

fulltax <- lookup %>%
  mutate(genome_id=as.character(genome_id)) %>%
  left_join(cleantax %>%
              rbind(blast_int)) %>%
  left_join(meta %>%
              select(seq_id,tiplab,length,wiki_link,both_ids))

kraken2 <- tbl(con,"kraken2_contigs") %>% collect()

#get metabolomics data
scfa <- read.csv(file="~/Documents/DFI/Projects/matthew/4_10_20_scfa.csv",
                 stringsAsFactors = F) %>%
  select(-X)

no_grow <-  c("MSK.11.31", "MSK.13.2", "MSK.14.9")

}



# load data ---------------------------------------------------------------


#save.image("/Volumes/pamer-lab/Eric.Littmann/shiny_apps/shiny_image2.RData")
#image #3 contains seq data.. 1.7Gb
#load("/Volumes/pamer-lab/Eric.Littmann/shiny_apps/shiny_image3.RData")
load("/pamer-lab/Eric.Littmann/shiny_apps/shiny_image4.RData")

#save.image("~/Desktop/shiny_apps/shiny_image4.RData")
#save.image("/Volumes/pamer-lab/Eric.Littmann/shiny_apps/shiny_image4.RData")

#system("ls -lah /Volumes/pamer-lab/Eric.Littmann/shiny_apps/shiny*3*")

scfa <- read.csv(file="/pamer-lab/Eric.Littmann/shiny_apps/4_10_20_scfa.csv",
                 stringsAsFactors = F) %>%
  select(-X)

no_grow <-  c("MSK.11.31", "MSK.13.2", "MSK.14.9")

#TM81 TM6 TM13

fulltax <- fulltax %>%
  unique()

# functions ---------------------------------------------------------------

make_rect_phylo <- function(bp_value=1300,
                            phyla_filter=c("Firmicutes","Actinobacteria","Bacteroidetes","Proteobacteria"),
                            show_species=F,
                            show_msk_id=F,
                            text_size=2.5){
  #bp_value=1300
  #phyla_filter=c("Firmicutes")
  #text_size=2.5
  to_remove2 <- meta %>%
    filter(phylum %!in% phyla_filter) %>%
    select(tiplab) %>%
    mutate(tiplab=as.character(tiplab))
  
  to_remove <- meta %>%
    filter(seq_id %!in% blast_top$seq_id|length < bp_value|seq_id %in% c("TM81","TM6","TM13",
                                                                          "TM31","TM32","MSK1_4","MSK1_28")) %>%
    select(tiplab) %>%
    mutate(tiplab=as.character(tiplab)) %>%
    rbind(to_remove2)
  
  physub <- drop.tip(tree,to_remove$tiplab)
  physub$edge.length[physub$edge.length<0] <- 0
  
  
  if(show_species==T & show_msk_id==F){
    gg <- ggtree(physub,
                 layout="rectangular",
                 ladderize=F,size=1) %<+%
      meta +
      scale_color_discrete(guide=F) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Bacteroidetes"),fill="skyblue",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Lachnospiraceae"),fill="#EC9B96",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Ruminococcaceae"),fill="#9AAE73",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Actinobacteria"),fill="#A77097",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Proteobacteria"),fill="red",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Erysipelotrichaceae"),fill="orange",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=class,value="Negativicutes"),fill="green4",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=class,value="Bacilli"),fill="blue",xend=1.5,xadd=1.5,size=6) +
      geom_tippoint(size=1,alpha=0.7) +
      #yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Tissierellaceae"),fill="gray4") +
      geom_tiplab(aes(label=species),
                  size=text_size,
                  hjust=-0.2,align=F) +
      theme(#legend.title=element_blank(),
        legend.position="top",
        plot.title=element_text(size=20,hjust = 0.5)) +
      xlim(c(0,4))
  }else if(show_msk_id==T & show_species==F){
    gg <- ggtree(physub,
                 layout="rectangular",
                 ladderize=F,size=1) %<+%
      meta +
      scale_color_discrete(guide=F) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Bacteroidetes"),fill="skyblue",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Lachnospiraceae"),fill="#EC9B96",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Ruminococcaceae"),fill="#9AAE73",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Actinobacteria"),fill="#A77097",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Proteobacteria"),fill="red",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Erysipelotrichaceae"),fill="orange",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=class,value="Negativicutes"),fill="green4",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=class,value="Bacilli"),fill="blue",xend=1.5,xadd=1.5,size=6) +
      geom_tippoint(size=1,alpha=0.7) +
      #yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Tissierellaceae"),fill="gray4") +
      geom_tiplab(aes(label=msk_id),
                  size=text_size,
                  hjust=-0.2,align=F) +
      theme(#legend.title=element_blank(),
        legend.position="top",
        plot.title=element_text(size=20,hjust = 0.5)) +
      xlim(c(0,4))
  }else if(show_msk_id==T & show_species==T){
    gg <- ggtree(physub,
                 layout="rectangular",
                 ladderize=F,size=1) %<+%
      meta +
      scale_color_discrete(guide=F) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Bacteroidetes"),fill="skyblue",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Lachnospiraceae"),fill="#EC9B96",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Ruminococcaceae"),fill="#9AAE73",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Actinobacteria"),fill="#A77097",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Proteobacteria"),fill="red",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Erysipelotrichaceae"),fill="orange",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=class,value="Negativicutes"),fill="green4",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=class,value="Bacilli"),fill="blue",xend=1.5,xadd=1.5,size=6) +
      geom_tippoint(size=1,alpha=0.7) +
      #yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Tissierellaceae"),fill="gray4") +
      geom_tiplab(aes(label=both_ids),
                  size=text_size,
                  hjust=-0.2,align=F) +
      theme(#legend.title=element_blank(),
        legend.position="top",
        plot.title=element_text(size=20,hjust = 0.5)) +
      xlim(c(0,4))
  }else{
    gg <- ggtree(physub,
                 layout="rectangular",
                 ladderize=F,size=1) %<+%
      meta +
      scale_color_discrete(guide=F) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Bacteroidetes"),fill="skyblue",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Lachnospiraceae"),fill="#EC9B96",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Ruminococcaceae"),fill="#9AAE73",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Actinobacteria"),fill="#A77097",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=phylum,value="Proteobacteria"),fill="red",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Erysipelotrichaceae"),fill="orange",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=class,value="Negativicutes"),fill="green4",xend=1.5,xadd=1.5,size=6) +
      yingtools2::geom_hilight(aes(isTip=isTip,var=class,value="Bacilli"),fill="blue",xend=1.5,xadd=1.5,size=6) +
      geom_tippoint(size=1,alpha=0.7) +
      #yingtools2::geom_hilight(aes(isTip=isTip,var=family,value="Tissierellaceae"),fill="gray4") +
      #geom_tiplab(aes(label=species),
      #            size=text_size,
      #            hjust=-0.2,align=F) +
      theme(#legend.title=element_blank(),
        legend.position="top",
        plot.title=element_text(size=20,hjust = 0.5)) +
      xlim(c(0,4))
  }
  #gg
  return(gg)
  
}

make_metabomoics <- function(show_msk_id=T){
  if(show_msk_id){
  gg <- scfa %>%
    inner_join(blast %>%
                 select(msk_id,seq_id,species)) %>%
    # inner_join(matt_tax %>%
    #             select(msk_id,seq_id,species=matt_species)) %>%
    reshape2::melt(id.vars=c("msk_id","seq_id","species")) %>%
    dplyr::rename(compound=variable) %>%
    ggplot(aes(x=msk_id,y=value,fill=compound)) +
    geom_bar(stat="identity") +
    # scale_fill_gradient2(trans=log_epsilon_trans(epsilon=1),
    #                      midpoint=5) +
    facet_grid(compound ~ species,scales="free",space="free_x") +
    theme(axis.text.x=element_text(angle=90),
          strip.text.x=element_text(angle=90),
          strip.text.y=element_text(angle=0),
          legend.position="none")
  }else{
    gg <- scfa %>%
      inner_join(blast %>%
                   select(msk_id,seq_id,species)) %>%
      reshape2::melt(id.vars=c("msk_id","seq_id","species")) %>%
      dplyr::rename(compound=variable) %>%
      ggplot(aes(x=msk_id,y=value,fill=compound)) +
      geom_bar(stat="identity") +
      # scale_fill_gradient2(trans=log_epsilon_trans(epsilon=1),
      #                      midpoint=5) +
      facet_grid(compound ~ species,scales="free",space="free_x") +
      theme(axis.text.x=element_blank(),
            strip.text.x=element_text(angle=90),
            strip.text.y=element_text(angle=0),
            legend.position="none")
  }
  return(gg)
  gg
}

make_kraken_plot <- function(mskid){
  
  int <- kraken2 %>%
    left_join(lookup %>%
                select(msk_id,seq_id)) %>%
    filter(msk_id==mskid) %>%
    group_by(seq_id) %>%
    summarize(total=sum(length)) 
  
  gg <- kraken2 %>%
    left_join(lookup %>%
                select(msk_id,seq_id)) %>%
    filter(msk_id==mskid) %>%
    group_by(seq_id) %>%
    mutate(total=sum(length)) %>%
    ungroup() %>%
    group_by(msk_id,seq_id,total,taxon) %>%
    summarize(length=sum(length)) %>%
    ungroup() %>%
    #left_join(fulltax) %>%
    mutate(pctseqs=length/total,
           # pctseqs=length,
           taxlabel=gsub("\\(.+","",taxon),
           taxlabel=gsub(" ","\n",taxlabel),
           tlabel=ifelse(pctseqs > 0.25,taxlabel,NA)) %>%
    group_by(seq_id) %>%
    arrange(taxon) %>%
    mutate(cum.pct=cumsum(pctseqs),
           prev=lag(cum.pct)) %>%
    replace_na(list(prev=0)) %>%
    mutate(y.text=(prev+cum.pct)/2) %>%
    ggplot() +
    geom_bar(aes(x=seq_id,y=pctseqs,fill=taxon),position="fill",stat="identity") +
    geom_text(aes(x=seq_id,y=1-y.text,label=tlabel)) +
    theme_bw() +
    theme(legend.position = "none") +
    #facet_wrap("species") +
    ylab("relative abundance") +
    ggtitle(paste("contig lengths per taxon added together\ntotal genome size:",int$total))
  
  return(gg)
  
  
}

# test area --------------------------------------------------------------------
if(F){

make_rect_phylo(bp_value=1300)
  
  make_kraken_plot("MSK.1.10")
  
#make a kraken2 function
glimpse(kraken2)

make_kraken_plot <- function(mskid){
  
  gg <- kraken2 %>%
    left_join(lookup %>%
                select(msk_id,seq_id)) %>%
    filter(msk_id==mskid) %>%
    group_by(seq_id) %>%
    mutate(total=sum(length)) %>%
    ungroup() %>%
    group_by(msk_id,seq_id,total,taxon) %>%
    summarize(length=sum(length)) %>%
    ungroup() %>%
    left_join(fulltax) %>%
    mutate(pctseqs=length/total,
           taxlabel=gsub("\\(.+","",taxon),
           taxlabel=gsub(" ","\n",taxlabel),
           tlabel=ifelse(pctseqs > 0.25,taxlabel,NA)) %>%
    group_by(seq_id) %>%
    arrange(taxon) %>%
    mutate(cum.pct=cumsum(pctseqs),
           prev=lag(cum.pct)) %>%
    replace_na(list(prev=0)) %>%
    mutate(y.text=(prev+cum.pct)/2) %>%
    ggplot() +
    geom_bar(aes(x=seq_id,y=pctseqs,fill=taxon),position="fill",stat="identity") +
    geom_text(aes(x=seq_id,y=y.text,label=tlabel)) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap("species")
  
  return(gg)
    
    
}

  
}
# ui ----------------------------------------------------------------------



ui <- fluidPage(
  theme = shinytheme("yeti"),
  titlePanel("DFI Commensal Library v4.6"),
  
    mainPanel(uiOutput("url"),
              DT::dataTableOutput("bug_list"),
              tabsetPanel(type="tabs",
                tabPanel("Phylogenetic Tree",
                         div(style="display: inline-block;vertical-align:top; width=100px;",
                             checkboxGroupInput("phyla","Phyla:",
                                                inline=T,
                                                c("Firmicutes","Actinobacteria","Bacteroidetes","Proteobacteria"),
                                                selected= c("Firmicutes","Actinobacteria","Bacteroidetes","Proteobacteria"))),
                         div(style="display: inline:block;vertical-align:top; width=100px;",
                             checkboxInput("msk_id","show msk_id",F)),
                         div(style="display: inline:block;vertical-align:top; width=100px;",
                             checkboxInput("species","show species",F)),
                         plotOutput("tree",height="1200px",width="1200px")),
                tabPanel("BLAST biobank protein",
                         selectInput("program", "Program:", choices=c("blastn","blastp","blastx"), width="100px",
                                     selected="blastp"),
                         textAreaInput('queryProt', 'Input sequence:', 
                                                       value = "", 
                                                       placeholder = "", 
                                                       width = "600px", height="200px"),
                         DT::dataTableOutput("blastProtResults"),
                         actionButton("blastProt", "Run BLAST"),
                         downloadButton("downloadProtData", "Download")),
                tabPanel("BLAST biobank 16s rRNA",
                         textAreaInput('query16s', 'Input 16s rRNA Sanger sequence:', 
                                       value = "", 
                                       placeholder = "", 
                                       width = "600px", height="200px"),
                         actionButton("blast16s", "Run BLASTn"),
                         DT::dataTableOutput("blast16sResults"),
                         downloadButton("download16sData", "Download")),
                tabPanel("kraken2 contig contam",
                         textInput("kraken_mskid","Enter msk_id:",value="MSK.1.10"),
                         plotOutput("kraken2")),
                tabPanel("Metabolomics",
                         checkboxInput("msk_id","show msk_id",F),
                         plotOutput("metabolomics",height="500px",width="1200px"))
             )
    )
    
)


# server ------------------------------------------------------------------
server <- function(input, output, session) {
  
  url <- a("DFI wiki",href="http://128.135.41.103/wiki/index.php")
  output$url <- renderUI({
    tagList("Wiki link:",url)
  })
  
  #tree plot
  tree_reactive <- reactive({
    make_rect_phylo(#bp_value=input$bp_cutoff,
                    bp_value=1300,
                    phyla_filter=input$phyla,
                    show_msk_id=input$msk_id,
                    show_species=input$species,
                    text_size = 2.5)
  })
  output$tree <- renderPlot({
    tree_reactive()
  })
  
  
  #BLAST protein panel ----------
 # db <- "/Volumes/pamer-lab/Eric.Littmann/Projects/sam_light/msk_library_prot/msk_library_all_prots"
  #db <- "/Volumes/pamer-lab/Eric.Littmann/Projects/sam_light/msk_library_prot/lux"
  protdb <- "/pamer-lab/Eric.Littmann/shiny_apps/blast_db/msk_library_all_prots"
  
  blastProtResults <- eventReactive(input$blastProt, {
    query <- input$queryProt
    
    tmp <- tempfile(fileext = ".fa")
    
    #make sure fasta is formatted properly
    if(startsWith(query,">")){
      writeLines(query, tmp)
    }else{
      writeLines(paste0(">Query\n",query), tmp)
    }
    
    print("running BLAST...")
    
    withProgress(message="running BLAST.. this may take a few minutes..",{
    bdata <- system(paste0(input$program,
                          #"blastp",
                          " -query ",
                          tmp," -db ",
                          protdb,
                          " -outfmt 6"),intern=T)
    })
    
    int <- bdata %>%
      as.data.frame() 
    colnames(int) <- "V1"
    int <- int %>%
      tidyr::separate(V1,into= c("qseqid","sseqid","pident",
                                 "length","mismatch",
                                 "gapopen","qstart",
                                 "qend","sstart","send",
                                 "evalue","bitscore"),
                      sep="\t")
    hits <- int %>%
      mutate(seq_id=gsub("\\_\\_.+","",sseqid),
             pident=as.numeric(pident),
             length=as.numeric(length),
             evalue=as.numeric(evalue),
             bitscore=as.numeric(bitscore)) %>%
      filter(pident >= 20) %>%
      left_join(lookup) %>%
      inner_join(blast %>%
                  select(seq_id,phylum,class,order,family,genus,species)) %>%
      ungroup() %>%
      group_by(seq_id,sseqid) %>%
      arrange(desc(bitscore)) %>%
      dplyr::slice(1) %>%
      #filter(!phylum=="Ascomycota") %>%
      separate(sseqid,into=c("seq_id","locus_tag"),sep="\\_\\_") %>%
      ungroup() %>%
      select(msk_id,seq_id,phylum,
             family,species,qseqid,
             locus_tag,pident,length,evalue,bitscore) %>%
      left_join(seq %>%
                  select(locus_tag,prot_sequence,nuc_sequence))
    
    if(nrow(hits)==0){
      hits <- c("no significant hits...") %>%
        as.data.frame()
    }else{
      hits <- hits
    }
  },ignoreNULL = T)
    

  output$blastProtResults <- DT::renderDataTable({
    if (is.null(blastProtResults())){
    } else {
      blastProtResults()
    }
  }, selection="single")
  
  output$downloadProtData <- downloadHandler(
    
    filename = function(){
      paste0("blast_results_",Sys.Date(),".csv")
    },
    
    content = function(file) {
      write.csv(blastProtResults(),file,
                  row.names=F,quote=F)
    }
  )

  # blast 16s rRNA gene ---------------
  db16s <- "/pamer-lab/EddiLin/database/biobank_16s/biobank.16s.all"
  
  blast16sResults <- eventReactive(input$blast16s, {
    query16s <- input$query16s
    
    tmp16s <- tempfile(fileext = ".fa")
    
    #make sure fasta is formatted properly
    if(startsWith(query16s,">")){
      writeLines(query16s, tmp16s)
    }else{
      writeLines(paste0(">Query\n",query16s), tmp16s)
    }
    
    print("running BLAST...")
    
    withProgress(message="running BLASTn.. this may take a few minutes..",{
      b16sdata <- system(paste0("blastn",
                             " -query ",
                             tmp16s," -db ",
                             db16s,
                             " -outfmt 6"),intern=T)
    })
    
    int16s <- b16sdata %>%
      as.data.frame() 
    colnames(int16s) <- "V1"
    int16s <- int16s %>%
      tidyr::separate(V1,into= c("qseqid","sseqid","pident",
                                 "length","mismatch",
                                 "gapopen","qstart",
                                 "qend","sstart","send",
                                 "evalue","bitscore"),
                      sep="\t", convert = T) %>%
      as_tibble()
    hits16s <- int16s %>%
      separate(sseqid, into=c("seq_id","locus_tag","lengthBP"), 
               sep = "-", remove = F) %>%
      filter(pident >= 50) %>%
      left_join(lookup) %>%
      inner_join(blast %>%
                   dplyr::select(seq_id,phylum,class,order,family,genus,species)) %>%
      ungroup() %>%
      group_by(sseqid) %>%
      top_n(1, wt = bitscore) %>%
      #filter(!phylum=="Ascomycota") %>%
      ungroup() %>%
      dplyr::select(msk_id,seq_id,phylum,
             family,species, qseqid,
             locus_tag, pident,length,evalue,bitscore) %>%
      left_join(seq %>%
                  dplyr::select(locus_tag,nuc_sequence))
    
    if(nrow(hits16s)==0){
      hits16s <- c("no significant hits...") %>%
        as.data.frame()
    }else{
      hits16s <- hits16s
    }
  },ignoreNULL = T)
  
  output$blast16sResults <- DT::renderDataTable({
    if (is.null(blast16sResults())){
    } else {
      blast16sResults()
    }
  }, selection="single")
  
  output$download16sData <- downloadHandler(
    
    filename = function(){
      paste0("blastn16s_results_",Sys.Date(),".csv")
    },
    
    content = function(file) {
      write.csv(blast16sResults(),file,
                row.names=F,quote=F)
    }
  )
  
  #metabolomics plot ----------------
  metabolomics_reactive <- reactive({
    make_metabomoics(show_msk_id=input$msk_id)
  })
  output$metabolomics <- renderPlot({
    metabolomics_reactive()
  })
  
  kraken2_reactive <- reactive({
    make_kraken_plot(mskid=input$kraken_mskid)
  })
  
  output$kraken2 <- renderPlot({
    kraken2_reactive()
  })
  
  #datatable ---------------
  # tbl <- meta %>%
  #   filter(!is.na(phylum)) %>%
  #   select(msk_id,seq_id,wiki_link) %>%
  #   left_join(blast %>%
  #               select(msk_id,seq_id,phylum,class,order,family,genus,species)) %>%
  #   filter(!phylum=="Ascomycota") %>%
  #   mutate(wiki_link=paste0("<a href=\"",wiki_link,
  #                           "\",target=\"_blank\">",
  #                           wiki_link,"</a>"))%>%
  #   select(msk_id,seq_id,phylum,class,order,family,
  #          genus,species,wiki_link) %>%
  #   arrange(msk_id)
  
  # tbl <- fulltax %>%
  #   select(msk_id,seq_id,phylum,class,order,family,
  #          genus,species,wiki_link) %>%
   
     tbl <- fulltax %>%
    select(msk_id,seq_id,phylum,class,order,family,
           genus,species,wiki_link) %>%
    mutate(wiki_link=paste0("<a href=\"",wiki_link,
                            "\",target=\"_blank\">",
                            wiki_link,"</a>")) %>%
    arrange(msk_id)
  
    #
  
  output$bug_list <- DT::renderDataTable(
    #datatable(tbl) %>%
      #formatStyle(colnames(tbl),background = "black"),
    tbl,
   options = list(pageLength=5),
   escape=F
  )
}

# Run app ----
#shinyApp(ui, server)
runApp(list(ui=ui,server=server),host="0.0.0.0",port=5050)
