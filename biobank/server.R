library(shiny)
library(tidyverse)
library(stringr)
library(yingtools2)
library(shinythemes)
library(ggtree)
library(grid)
library(gridExtra)
library(DT)
library(ape)

# Sys.setenv(SHELL = "/bin/bash")

# load data ---------------------------------------------------------------

# change this flag to FALSE when running on virtual box
localF = T

if (localF){
  # local 
  load("~/OneDrive - The University of Chicago/dfimmfshiny/biobank/shiny_image_sm.RData")
  
  scfa <- read.csv(file="/Volumes/pamer-lab/Eric.Littmann/shiny_apps/4_10_20_scfa.csv",
                   stringsAsFactors = F) %>%
    select(-X)
  
  protdb <- "/Volumes/pamer-lab/Eric.Littmann/shiny_apps/blast_db/msk_library_all_prots"
  
  db16s <- "/Volumes/pamer-lab/EddiLin/database/biobank_16s/biobank.16s.all"

} else {
  # virtual machine version
   load("/srv/shiny-server/biobank/shiny_image_sm.RData")
#  load("/pamer-lab/Eric.Littmann/shiny_apps/shiny_image_sm.RData") 
#  load("/pamer-lab/Eric.Littmann/shiny_apps/shiny_image4.RData")

scfa <- read.csv(file="/srv/shiny-server/biobank/4_10_20_scfa.csv",
                   stringsAsFactors = F) %>%
    select(-X) 

 
#  scfa <- read.csv(file="/pamer-lab/Eric.Littmann/shiny_apps/4_10_20_scfa.csv",
#                   stringsAsFactors = F) %>%
#    select(-X)
  
#  protdb <- "/pamer-lab/Eric.Littmann/shiny_apps/blast_db/msk_library_all_prots"
protdb <- "/srv/shiny-server/biobank/msk_library_all_prots"  
#  db16s <- "/pamer-lab/EddiLin/database/biobank_16s/biobank.16s.all"
db16s <- "/srv/shiny-server/biobank/biobank.16s.all"

  
}


#### END of switch ####
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
  
  intk <- kraken2 %>%
    left_join(lookup %>%
                select(msk_id,seq_id)) %>%
    filter(msk_id %in% mskid) %>%
    group_by(seq_id,msk_id) %>%
    summarize(total=sum(length)) %>%
    mutate(sideLab = paste0("MSK_id:",msk_id,"\nTotal bp:", total))
  
  bardf <- kraken2 %>%
    left_join(lookup %>%
                select(msk_id,seq_id)) %>%
    filter(msk_id %in% mskid) %>%
    group_by(msk_id, seq_id) %>%
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
    mutate(y.text=(prev+cum.pct)/2) 
  
  gg <- ggplot(bardf) +
    geom_bar(aes(x=seq_id,y=pctseqs,fill=taxon),
             position="fill",stat="identity",
             width = 0.75) +
    geom_text(aes(x=seq_id,y=1-y.text,label=tlabel)) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "contig lengths per taxon added together",
         y = "relative abundance") +
    geom_text(data = intk, 
              aes(x= seq_id, y = 0.5, label = sideLab),
              hjust = 0.5, nudge_x = -0.45, angle = 90) 
  
  return(gg)
  
  
}


# server ------------------------------------------------------------------

shinyServer(function(input, output, session) {
  
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
      bdata <- system(paste0("/usr/local/sbin/",input$program,
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
      b16sdata <- system(paste0("/usr/local/sbin/blastn",
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
    make_metabomoics(show_msk_id=input$metamsk_id)
  })
  output$metabolomics <- renderPlot({
    metabolomics_reactive()
  })

  
  kraken2_reactive <- reactive({
    kids <- trimws(str_split(input$kraken_mskid,",")[[1]])
    make_kraken_plot(mskid = kids)
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
    class = "display nowrap compact",
    options = list(pageLength=25, scrollX = TRUE,
                   search = list(regex = T, caseInsensitive = T)),
    escape = F,
    filter = "top"
  )
})
